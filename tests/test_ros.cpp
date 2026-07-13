// ROS bridge: the node graph delivers published messages to every subscriber, and
// the sim→ROS adapters produce well-formed sensor_msgs / tf2 messages.
#include "phys/ros_bridge.h"
#include "check.h"
#include <cmath>
using namespace phys::ros;

int main() {
    Graph::get().reset();
    Node node("physics_sim");

    // A) Imu publish → subscribe delivers the exact message
    {
        auto pub = node.create_publisher<Imu>("/imu");
        Imu got; int calls = 0;
        auto sub = node.create_subscription<Imu>("/imu", [&](const Imu& m) { got = m; calls++; });
        Imu msg; msg.header.frame_id = "imu_link"; msg.linear_acceleration = {0, 9.81, 0}; msg.angular_velocity = {0.1, 0, -0.2};
        pub.publish(msg);
        CHECK(calls == 1);
        CHECK_NEAR(got.linear_acceleration.y, 9.81, 1e-9);
        CHECK_NEAR(got.angular_velocity.z, -0.2, 1e-9);
        CHECK(got.header.frame_id == "imu_link");
    }

    // B) two subscribers both receive a LaserScan; ranges round-trip
    {
        auto pub = node.create_publisher<LaserScan>("/scan");
        int a = 0, b = 0; float lastRange = 0;
        auto s1 = node.create_subscription<LaserScan>("/scan", [&](const LaserScan& m) { a++; lastRange = m.ranges.empty() ? 0 : m.ranges[1]; });
        auto s2 = node.create_subscription<LaserScan>("/scan", [&](const LaserScan&) { b++; });
        LaserScan scan; scan.ranges = {1.0f, 2.5f, 3.0f}; scan.range_max = 5.0f;
        pub.publish(scan); pub.publish(scan);
        CHECK(a == 2 && b == 2);
        CHECK_NEAR(lastRange, 2.5, 1e-6);
        CHECK(Graph::get().pubCount["/scan"] == 2);
    }

    // C) JointState round-trips names + positions
    {
        auto pub = node.create_publisher<JointState>("/joint_states");
        JointState got;
        auto sub = node.create_subscription<JointState>("/joint_states", [&](const JointState& m) { got = m; });
        JointState js; js.name = {"j0", "j1"}; js.position = {0.5, -1.2}; js.velocity = {0.0, 0.3};
        pub.publish(js);
        CHECK(got.name.size() == 2 && got.name[1] == "j1");
        CHECK_NEAR(got.position[0], 0.5, 1e-9);
        CHECK_NEAR(got.velocity[1], 0.3, 1e-9);
    }

    // D) a publisher with no subscribers is a no-op (doesn't crash), count still increments
    {
        auto pub = node.create_publisher<Imu>("/unheard");
        pub.publish(Imu{});
        CHECK(Graph::get().pubCount["/unheard"] == 1);
    }

    // E) sim → ROS adapters build correct messages
    {
        Imu m = makeImu(phys::Vector3(0, 9.81, 0), phys::Vector3(0, 0, 0.4));
        CHECK_NEAR(m.linear_acceleration.y, 9.81, 1e-9);
        CHECK_NEAR(m.angular_velocity.z, 0.4, 1e-9);

        LaserScan s = makeScan({2.0f, 2.0f, 2.0f, 2.0f, 2.0f}, 3.14159f, 6.0f);
        CHECK(s.ranges.size() == 5);
        CHECK_NEAR(s.angle_min, -3.14159 / 2, 1e-4);
        CHECK_NEAR(s.angle_increment, 3.14159 / 4, 1e-4);

        phys::Quaternion q(1, 0, 0, 0);
        TransformStamped t = makeTf(phys::Vector3(1, 2, 3), q, "world", "base");
        CHECK(t.child_frame_id == "base" && t.header.frame_id == "world");
        CHECK_NEAR(t.transform.translation.x, 1.0, 1e-9);
        CHECK_NEAR(t.transform.rotation.w, 1.0, 1e-9);            // ROS order is (x,y,z,w)
    }

    return test::report("ros");
}
