// Robotics layer: URDF/MJCF loading, IMU / Lidar / contact-force sensors,
// recursive Newton–Euler inverse dynamics, Jacobian-transpose IK and tendon /
// muscle actuators — each checked against an analytic or known-physical result.
#include "phys/loader.h"
#include "phys/robotics.h"
#include "check.h"
#include <cmath>
using namespace phys;

int main() {
    // A) URDF loader — 2 links + 1 revolute joint with the right masses & axis.
    {
        std::string urdf = R"URDF(
        <robot name="arm">
          <link name="base">
            <inertial>
              <mass value="1.0"/>
              <inertia ixx="0.10" iyy="0.10" izz="0.05"/>
            </inertial>
            <collision><geometry><cylinder radius="0.05" length="0.2"/></geometry></collision>
          </link>
          <link name="link1">
            <inertial>
              <mass value="2.0"/>
              <inertia ixx="0.20" iyy="0.20" izz="0.02"/>
            </inertial>
            <collision><geometry><box size="0.1 0.1 0.6"/></geometry></collision>
          </link>
          <joint name="shoulder" type="revolute">
            <parent link="base"/>
            <child link="link1"/>
            <axis xyz="0 0 1"/>
            <origin xyz="0 0 0.3"/>
          </joint>
        </robot>)URDF";
        Model m = loadURDF(urdf);
        CHECK(m.links.size() == 2);
        CHECK(m.joints.size() == 1);
        CHECK(m.links[0].name == "base");
        CHECK(m.links[1].name == "link1");
        CHECK_NEAR(m.links[0].mass, 1.0, 1e-9);
        CHECK_NEAR(m.links[1].mass, 2.0, 1e-9);
        CHECK(m.links[0].geom == GeomType::Cylinder);
        CHECK(m.links[1].geom == GeomType::Box);
        CHECK_NEAR(m.links[1].size.z, 0.6, 1e-9);            // box extents parsed
        CHECK(m.joints[0].type == JointKind::Revolute);
        CHECK(m.joints[0].parent == "base");
        CHECK(m.joints[0].child == "link1");
        CHECK_NEAR(m.joints[0].axis.x, 0.0, 1e-9);
        CHECK_NEAR(m.joints[0].axis.z, 1.0, 1e-9);           // axis "0 0 1"
        CHECK_NEAR(m.joints[0].origin.z, 0.3, 1e-9);
        CHECK(m.findLink("link1") != nullptr);
    }

    // A2) MJCF loader — nesting recovers the parent/child tree.
    {
        std::string mjcf = R"MJCF(
        <mujoco>
          <worldbody>
            <body name="upper" pos="0 0 0.5">
              <joint name="j0" type="hinge" axis="0 1 0"/>
              <geom type="box" size="0.05 0.05 0.25" mass="1.5"/>
              <body name="fore" pos="0 0 0.5">
                <joint name="j1" type="slide" axis="0 0 1"/>
                <geom type="sphere" size="0.04" mass="0.7"/>
              </body>
            </body>
          </worldbody>
        </mujoco>)MJCF";
        Model m = loadMJCF(mjcf);
        CHECK(m.links.size() == 2);
        CHECK(m.joints.size() == 2);
        CHECK_NEAR(m.links[0].mass, 1.5, 1e-9);
        CHECK_NEAR(m.links[1].mass, 0.7, 1e-9);
        CHECK(m.joints[0].parent == "world");
        CHECK(m.joints[1].parent == "upper");               // nesting → parent
        CHECK(m.joints[1].child == "fore");
        CHECK(m.joints[0].type == JointKind::Revolute);      // hinge
        CHECK(m.joints[1].type == JointKind::Prismatic);     // slide
        CHECK_NEAR(m.joints[0].axis.y, 1.0, 1e-9);
    }

    // B) IMU on a static body: |accel| ≈ g pointing up, zero gyro.
    {
        IMU imu;                                             // gravity (0,-9.81,0)
        IMU::Reading r = imu.measure(Vector3(), Vector3(), Quaternion());  // at rest, identity
        CHECK_NEAR(r.accel.magnitude(), 9.81, 1e-9);
        CHECK_NEAR(r.accel.x, 0.0, 1e-9);
        CHECK_NEAR(r.accel.y, 9.81, 1e-9);                   // upward (reaction to gravity)
        CHECK_NEAR(r.accel.z, 0.0, 1e-9);
        CHECK_NEAR(r.gyro.magnitude(), 0.0, 1e-12);
        // rotated 90° about x: "up" now reads along the sensor's rotated axis, |·| still g.
        Quaternion tilt(std::cos(real_pi / 4), std::sin(real_pi / 4), 0, 0);
        IMU::Reading rt = imu.measure(Vector3(), Vector3(0, 0, 0), tilt);
        CHECK_NEAR(rt.accel.magnitude(), 9.81, 1e-9);
    }

    // C) RNE: a 1-link horizontal arm needs τ = m·g·l for the gravity term.
    {
        real m = 2.0, g = 9.81, l = 0.5;
        SerialChain c; c.gravity = Vector3(0, -g, 0);
        ChainLink L; L.axis = Vector3(0, 0, 1); L.com = Vector3(l, 0, 0);
        L.toChild = Vector3(2 * l, 0, 0); L.mass = m;
        c.links.push_back(L);
        std::vector<real> tau = c.inverseDynamics({0.0}, {0.0}, {0.0});
        CHECK(tau.size() == 1);
        CHECK_NEAR(tau[0], m * g * l, 1e-9);                 // gravity torque
        // Vertical (hanging straight down) → zero gravity torque about the axis.
        std::vector<real> tauV = c.inverseDynamics({-real_pi / 2}, {0.0}, {0.0});
        CHECK_NEAR(tauV[0], 0.0, 1e-9);
    }

    // D) Jacobian-transpose IK on a 2-link planar arm: error decreases each step
    //    (monotone) and reaches near zero.
    {
        SerialChain c;
        ChainLink a; a.axis = Vector3(0, 0, 1); a.toChild = Vector3(1, 0, 0); a.com = Vector3(0.5, 0, 0);
        c.links.push_back(a);
        c.links.push_back(a);                                // two unit links
        JacobianIK ik; ik.chain = &c;
        std::vector<real> q = {0.2, 0.2};
        Vector3 target(1.0, 0.6, 0.0);                       // |target| < 2 ⇒ reachable
        std::vector<real> hist;
        real finalErr = ik.solve(q, target, 400, &hist);
        bool monotone = true;
        for (size_t i = 1; i < hist.size(); i++)
            if (hist[i] > hist[i - 1] + 1e-12) monotone = false;
        CHECK(monotone);                                     // never got worse
        CHECK(hist.front() > 0.1);                           // it started far away
        CHECK(finalErr < 1e-3);                              // and converged
        CHECK_NEAR((c.endEffector(q) - target).magnitude(), 0.0, 1e-3);
    }

    // E) Lidar against a mock wall perpendicular to `forward` at distance d.
    {
        real d = 3.0;
        Lidar lid; lid.numRays = 5; lid.fov = real_pi / 2;   // ±45°, 5 rays
        lid.origin = Vector3(0, 0, 0);
        lid.forward = Vector3(1, 0, 0); lid.up = Vector3(0, 0, 1);   // fan in the x-y plane
        auto ray = [d](Vector3 o, Vector3 dir) -> real {
            if (dir.x <= 1e-9) return 100.0;                 // ray points away from the wall
            real t = (d - o.x) / dir.x;                      // hit the plane x = d
            return t > 0 ? t : 100.0;
        };
        std::vector<real> ranges = lid.scan(ray);
        CHECK(ranges.size() == 5);
        CHECK_NEAR(ranges[2], d, 1e-9);                      // centre ray → straight to the wall
        real edge = d / std::cos(real_pi / 4);               // outer rays at ±45°
        CHECK_NEAR(ranges[0], edge, 1e-6);
        CHECK_NEAR(ranges[4], edge, 1e-6);
        CHECK(ranges[1] > d && ranges[1] < edge);            // intermediate ray in between
    }

    // F) ContactForceSensor: accumulated impulse / dt = force.
    {
        ContactForceSensor s;
        s.addImpulse(Vector3(0, 10, 0));
        s.addImpulse(Vector3(2, 0, 0));
        Vector3 f = s.force(0.5);                            // (2,10,0) over 0.5 s
        CHECK_NEAR(f.x, 4.0, 1e-9);
        CHECK_NEAR(f.y, 20.0, 1e-9);
        s.clear();
        CHECK_NEAR(s.force(0.5).magnitude(), 0.0, 1e-12);
    }

    // G) Tendon (spring-damper cable) and Hill muscle actuators.
    {
        Tendon t; t.restLength = 1.0; t.stiffness = 100; t.damping = 0;
        CHECK_NEAR(t.tension(1.5, 0.0), 50.0, 1e-9);         // stretched 0.5 m ⇒ 50 N
        CHECK_NEAR(t.tension(0.5, 0.0), 0.0, 1e-9);          // slack ⇒ pulls nothing (no push)
        // routed length from waypoints
        Tendon r; r.route = {Vector3(0, 0, 0), Vector3(3, 4, 0), Vector3(3, 4, 12)};
        CHECK_NEAR(r.length(), 5.0 + 12.0, 1e-9);

        HillMuscle mus; mus.Fmax = 1000; mus.optLength = 0.1; mus.activation = 1;
        CHECK_NEAR(mus.tension(0.1, 0), 1000.0, 1e-6);       // peak active force at optimal length
        CHECK(mus.tension(0.1, 0) > mus.tension(0.13, 0));   // force-length bell falls off the peak
        mus.activation = 0.5;
        CHECK_NEAR(mus.tension(0.1, 0), 500.0, 1e-6);        // scales with neural drive
    }

    return test::report("robotics");
}
