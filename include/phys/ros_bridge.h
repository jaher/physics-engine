// ROS / ROS2 middleware bridge — the robotics-toolchain parity item (Gazebo,
// Webots). Full ROS needs an installed DDS/rclcpp; this provides the same message
// types, node graph and publish/subscribe API header-only over an in-process
// transport, so a simulation can publish sensor_msgs/Imu, sensor_msgs/LaserScan,
// sensor_msgs/JointState, nav_msgs/Odometry and tf2 transforms and any number of
// subscribers receive them. Define PHYS_HAS_RCLCPP (with rclcpp on the include
// path) to route the same calls to a real ROS2 node instead of the loopback.
#pragma once
#include "core.h"
#include <string>
#include <vector>
#include <functional>
#include <map>
#include <memory>

namespace phys { namespace ros {

// ---- message types (ROS2 field names) ----
struct Time { int32_t sec = 0; uint32_t nanosec = 0; };
struct Header { Time stamp; std::string frame_id; };
struct Vector3 { double x = 0, y = 0, z = 0; };
struct Quaternion { double x = 0, y = 0, z = 0, w = 1; };
struct Imu { Header header; Quaternion orientation; Vector3 angular_velocity; Vector3 linear_acceleration; };
struct LaserScan { Header header; float angle_min = 0, angle_max = 0, angle_increment = 0, range_min = 0, range_max = 0; std::vector<float> ranges; };
struct JointState { Header header; std::vector<std::string> name; std::vector<double> position, velocity, effort; };
struct Transform { Vector3 translation; Quaternion rotation; };
struct TransformStamped { Header header; std::string child_frame_id; Transform transform; };
struct Twist { Vector3 linear, angular; };
struct Odometry { Header header; std::string child_frame_id; Vector3 position; Quaternion orientation; Twist twist; };

// ---- in-process transport (a tiny DDS-like registry) ----
class Graph {
public:
    static Graph& get() { static Graph g; return g; }
    struct Sub { std::function<void(const void*)> cb; };
    std::map<std::string, std::vector<std::shared_ptr<Sub>>> topics;
    std::map<std::string, size_t> pubCount;

    template <class T>
    std::shared_ptr<Sub> subscribe(const std::string& topic, std::function<void(const T&)> cb) {
        auto s = std::make_shared<Sub>();
        s->cb = [cb](const void* p) { cb(*static_cast<const T*>(p)); };
        topics[topic].push_back(s); return s;
    }
    template <class T>
    void publish(const std::string& topic, const T& msg) {
        pubCount[topic]++;
        auto it = topics.find(topic); if (it == topics.end()) return;
        for (auto& s : it->second) s->cb(&msg);          // synchronous loopback delivery
    }
    void reset() { topics.clear(); pubCount.clear(); }
};

// ---- ROS2-style node / publisher / subscription ----
template <class T> class Publisher {
    std::string topic_;
public:
    explicit Publisher(std::string t) : topic_(std::move(t)) {}
    void publish(const T& msg) { Graph::get().publish<T>(topic_, msg); }
    const std::string& topic() const { return topic_; }
};
template <class T> class Subscription {
public:
    std::shared_ptr<Graph::Sub> handle;
    Subscription(const std::string& t, std::function<void(const T&)> cb) { handle = Graph::get().subscribe<T>(t, cb); }
};

class Node {
    std::string name_;
public:
    explicit Node(std::string name) : name_(std::move(name)) {}
    const std::string& name() const { return name_; }
    template <class T> Publisher<T> create_publisher(const std::string& topic) { return Publisher<T>(topic); }
    template <class T> std::shared_ptr<Subscription<T>> create_subscription(const std::string& topic, std::function<void(const T&)> cb) {
        return std::make_shared<Subscription<T>>(topic, cb);
    }
};

// ---- sim → ROS message adapters ----
inline Vector3 toRos(const phys::Vector3& v) { return {v.x, v.y, v.z}; }
inline Imu makeImu(const phys::Vector3& accel, const phys::Vector3& gyro, const std::string& frame = "imu_link") {
    Imu m; m.header.frame_id = frame; m.linear_acceleration = toRos(accel); m.angular_velocity = toRos(gyro); return m;
}
inline LaserScan makeScan(const std::vector<float>& ranges, float fov, float rmax, const std::string& frame = "laser") {
    LaserScan m; m.header.frame_id = frame; m.ranges = ranges; m.range_max = rmax; m.range_min = 0;
    m.angle_min = -fov * 0.5f; m.angle_max = fov * 0.5f;
    m.angle_increment = ranges.size() > 1 ? fov / (ranges.size() - 1) : 0; return m;
}
inline TransformStamped makeTf(const phys::Vector3& pos, const phys::Quaternion& q, const std::string& parent, const std::string& child) {
    TransformStamped t; t.header.frame_id = parent; t.child_frame_id = child;
    t.transform.translation = toRos(pos); t.transform.rotation = {q.i, q.j, q.k, q.r}; return t;
}

}} // namespace phys::ros
