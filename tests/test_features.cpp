// Tests for the framework-parity features: capsules, raycasts, joints + motors,
// heightfield, CCD, events, triggers, kinematic bodies, RK4, vehicle, character.
#include "phys/phys.h"
#include "check.h"
using namespace phys;

static RigidBody* mk(RigidBody& b, Vector3 pos, real mass = 1) {
    b.setMass(mass); Matrix3 I; I.setBlockInertiaTensor(Vector3(0.5, 0.5, 0.5), mass); b.setInertiaTensor(I);
    b.setPosition(pos); b.setDamping(0.9, 0.9); b.setCanSleep(false); b.calculateDerivedData(); return &b;
}

int main() {
    Contact contacts[8]; CollisionData data; data.contactArray = contacts;

    // ---- capsule collisions ----
    {
        RigidBody rb; mk(rb, Vector3(0, 0.75, 0));
        CollisionCapsule cap; cap.body = &rb; cap.radius = 0.3; cap.halfHeight = 0.5; cap.calculateInternals();
        CollisionPlane ground; ground.direction = Vector3(0, 1, 0); ground.offset = 0;
        data.reset(8); CapsuleDetector::capsuleAndHalfSpace(cap, ground, &data);
        CHECK(data.contactCount == 1);                      // only the lower end-sphere touches
        CHECK_NEAR(contacts[0].penetration, 0.05, 1e-9);    // 0.75-0.5-0.3
        RigidBody rs; mk(rs, Vector3(0.5, 0.8, 0));
        CollisionSphere sph; sph.body = &rs; sph.radius = 0.3; sph.calculateInternals();
        data.reset(8); CHECK(CapsuleDetector::capsuleAndSphere(cap, sph, &data) == 1);
        CHECK_NEAR(contacts[0].penetration, 0.1, 1e-9);     // 0.3+0.3-0.5
        RigidBody rc2; mk(rc2, Vector3(0.55, 0.8, 0));
        CollisionCapsule cap2; cap2.body = &rc2; cap2.radius = 0.3; cap2.halfHeight = 0.5; cap2.calculateInternals();
        data.reset(8); CHECK(CapsuleDetector::capsuleAndCapsule(cap, cap2, &data) == 1);
        CHECK_NEAR(contacts[0].penetration, 0.05, 1e-9);
        RigidBody rbx; mk(rbx, Vector3(0.65, 0.8, 0));
        CollisionBox box; box.body = &rbx; box.halfSize = Vector3(0.4, 0.4, 0.4); box.calculateInternals();
        data.reset(8); CHECK(CapsuleDetector::capsuleAndBox(cap, box, &data) >= 1);
    }

    // ---- raycasts ----
    {
        Ray r{Vector3(-5, 1, 0), Vector3(1, 0, 0)}; RayHit h;
        CHECK(raySphere(r, Vector3(0, 1, 0), 0.5, h)); CHECK_NEAR(h.t, 4.5, 1e-9);
        CHECK(rayPlane(Ray{Vector3(0, 5, 0), Vector3(0, -1, 0)}, Vector3(0, 1, 0), 0, h)); CHECK_NEAR(h.t, 5.0, 1e-9);
        RigidBody rb; mk(rb, Vector3(0, 1, 0));
        CollisionBox box; box.body = &rb; box.halfSize = Vector3(1, 1, 1); box.calculateInternals();
        CHECK(rayBox(Ray{Vector3(-5, 1, 0), Vector3(1, 0, 0)}, box, h)); CHECK_NEAR(h.t, 4.0, 1e-9);
        CHECK_NEAR(h.normal.x, -1.0, 1e-9);
        CollisionCapsule cap; cap.body = &rb; cap.radius = 0.4; cap.halfHeight = 0.5; cap.calculateInternals();
        CHECK(rayCapsule(Ray{Vector3(-5, 1, 0), Vector3(1, 0, 0)}, cap, h)); CHECK_NEAR(h.t, 4.6, 1e-3);
    }

    // ---- hinge: axis alignment + limits + motor ----
    {
        RigidBody a, b; mk(a, Vector3(0, 0, 0)); a.setInverseMass(0); a.setInverseInertiaTensor(Matrix3()); a.calculateDerivedData();
        mk(b, Vector3(1, 0, 0));
        HingeJoint hinge; hinge.set(&a, Vector3(0.5, 0, 0), Vector3(0, 0, 1), &b, Vector3(-0.5, 0, 0), Vector3(0, 0, 1));
        hinge.hasMotor = true; hinge.motorTargetVel = -2.0;   // spin b about z
        JointSolver solver; solver.joints.push_back(&hinge);
        for (int i = 0; i < 400; i++) { b.integrate(0.005); solver.solve(0.005); }
        CHECK(real_abs(hinge.angle()) > 0.5);                 // motor turned the joint
        Vector3 anchorGap = a.getPointInWorldSpace(Vector3(0.5, 0, 0)) - b.getPointInWorldSpace(Vector3(-0.5, 0, 0));
        CHECK(anchorGap.magnitude() < 0.02);                  // pivot held
        // now clamp with limits
        hinge.hasLimits = true; hinge.minAngle = -0.5; hinge.maxAngle = 0.5;
        for (int i = 0; i < 400; i++) { b.integrate(0.005); solver.solve(0.005); }
        CHECK(hinge.angle() > -0.62 && hinge.angle() < 0.62); // held at/near the limit
    }

    // ---- slider: off-axis lock + motor ----
    {
        RigidBody a, b; mk(a, Vector3(0, 0, 0)); a.setInverseMass(0); a.setInverseInertiaTensor(Matrix3()); a.calculateDerivedData();
        mk(b, Vector3(0, 0, 1));
        SliderJoint slider; slider.set(&a, Vector3(0, 0, 0), Vector3(0, 0, 1), &b, Vector3(0, 0, 0), Vector3(0, 0, 1));
        slider.hasMotor = true; slider.motorTargetVel = 1.0;
        JointSolver solver; solver.joints.push_back(&slider);
        b.setAcceleration(0, -9.81, 0);                       // gravity pulls off-axis
        for (int i = 0; i < 300; i++) { b.integrate(0.005); solver.solve(0.005); }
        CHECK(real_abs(b.getPosition().x) < 0.02);            // stayed on the rail
        CHECK(real_abs(b.getPosition().y) < 0.05);
        CHECK(b.getPosition().z > 1.5);                       // motor drove it along
    }

    // ---- distance joint holds length ----
    {
        RigidBody a, b; mk(a, Vector3(0, 0, 0)); a.setInverseMass(0); a.setInverseInertiaTensor(Matrix3()); a.calculateDerivedData();
        mk(b, Vector3(2, 0, 0)); b.setAcceleration(0, -9.81, 0);
        DistanceJoint dj; dj.set(&a, Vector3(), &b, Vector3(), 2.0);
        JointSolver solver; solver.joints.push_back(&dj);
        real maxErr = 0;
        for (int i = 0; i < 600; i++) { b.integrate(0.005); solver.solve(0.005);
            maxErr = std::max(maxErr, real_abs((b.getPosition() - a.getPosition()).magnitude() - 2.0)); }
        CHECK(maxErr < 0.05);
    }

    // ---- fixed joint locks pose ----
    {
        RigidBody a, b; mk(a, Vector3(0, 0, 0)); a.setInverseMass(0); a.setInverseInertiaTensor(Matrix3()); a.calculateDerivedData();
        mk(b, Vector3(1, 0, 0)); b.setAcceleration(0, -9.81, 0); b.setVelocity(0, 0, 3); b.setRotation(2, 1, 0.5);
        FixedJoint fj; fj.set(&a, Vector3(1, 0, 0), &b, Vector3(0, 0, 0));
        JointSolver solver; solver.iterations = 16; solver.joints.push_back(&fj);
        for (int i = 0; i < 500; i++) { b.integrate(0.005); solver.solve(0.005); }
        CHECK((b.getPosition() - Vector3(1, 0, 0)).magnitude() < 0.05);
        Matrix3 R; b.getOrientation(&R);
        CHECK_NEAR(R.data[0], 1.0, 0.05);                    // orientation stayed ≈ identity
    }

    // ---- gear couples angular velocities ----
    {
        RigidBody a, b; mk(a, Vector3(0, 0, 0)); mk(b, Vector3(2, 0, 0));
        a.setRotation(0, 0, 4); b.setRotation(0, 0, 0);
        GearJoint gear; gear.set(&a, Vector3(0, 0, 1), &b, Vector3(0, 0, 1), 2.0);
        JointSolver solver; solver.joints.push_back(&gear);
        for (int i = 0; i < 200; i++) { a.integrate(0.005); b.integrate(0.005); solver.solve(0.005); }
        real w0 = a.getRotation().z, w1 = b.getRotation().z;
        CHECK_NEAR(w0 + 2.0 * w1, 0.0, 0.02);                // constraint satisfied
        CHECK(real_abs(w1) > 0.1);                            // b actually driven
    }

    // ---- heightfield ----
    {
        HeightField hf; hf.init(11, 11, 1.0, Vector3(-5, 0, -5));
        for (int z = 0; z < 11; z++) for (int x = 0; x < 11; x++) hf.at(x, z) = 0.2 * x;   // ramp in x
        CHECK_NEAR(hf.sample(0, 0), 1.0, 1e-6);              // x=0 → sample idx 5 → 1.0
        Vector3 n = hf.normal(0, 0);
        CHECK(n.y > 0.9 && n.x < -0.1);                       // sloping up in +x → normal tilts -x
        RigidBody rb; mk(rb, Vector3(0, 1.2, 0));
        CollisionSphere s; s.body = &rb; s.radius = 0.5; s.calculateInternals();
        Contact cc[4]; CollisionData d2; d2.contactArray = cc; d2.reset(4);
        CHECK(hf.sphereContact(s, &d2) == 1);                 // 1.2 - 1.0 < 0.5 → contact
        Vector3 hit;
        CHECK(hf.raycast(Vector3(0, 5, 0), Vector3(0, -1, 0), 10, hit));
        CHECK_NEAR(hit.y, 1.0, 1e-3);
    }

    // ---- CCD: a bullet-fast sphere cannot tunnel a thin wall plane ----
    {
        RigidBody b; mk(b, Vector3(0, 1, -10), 1); b.setDamping(1, 1);
        b.setVelocity(0, 0, 400);                             // 400 m/s → 3.3 m per 1/120 step
        CCDWorld ccd; ccd.planes.push_back({Vector3(0, 0, -1), 0});   // wall at z=0 facing -z
        bool hit = false;
        for (int i = 0; i < 20 && !hit; i++) hit = ccd.integrate(b, 0.2, 1.0 / 120, 0.5);
        CHECK(hit);
        CHECK(b.getPosition().z <= 0.01);                     // stopped at the wall, no tunnel
        CHECK(b.getVelocity().z < 0);                          // reflected
    }

    // ---- contact events + trigger volumes ----
    {
        RigidBody a, b; mk(a, Vector3(0, 0, 0)); mk(b, Vector3(5, 0, 0));
        Contact one; one.body[0] = &a; one.body[1] = &b;
        int begins = 0, ends = 0;
        ContactEvents ev;
        ev.onBegin = [&](RigidBody*, RigidBody*) { begins++; };
        ev.onEnd = [&](RigidBody*, RigidBody*) { ends++; };
        ev.frame(&one, 1); ev.frame(&one, 1); ev.frame(nullptr, 0);
        CHECK(begins == 1); CHECK(ends == 1);
        TriggerVolume tv; tv.centre = Vector3(0, 0, 0); tv.isSphere = true; tv.radius = 1;
        int enters = 0, exits = 0;
        tv.onEnter = [&](RigidBody*) { enters++; }; tv.onExit = [&](RigidBody*) { exits++; };
        std::vector<RigidBody*> bodies{&a, &b};
        tv.frame(bodies);                                     // a inside
        a.setPosition(10, 0, 0); tv.frame(bodies);            // a left
        CHECK(enters == 1); CHECK(exits == 1);
    }

    // ---- kinematic body pushes a dynamic one (infinite mass + velocity) ----
    {
        RigidBody kin; mk(kin, Vector3(-2, 0.5, 0)); makeKinematic(kin);
        CHECK(kin.getInverseMass() == 0);
        moveKinematic(kin, Vector3(2, 0, 0), 0.5);
        CHECK_NEAR(kin.getPosition().x, -1.0, 1e-9);
        CHECK_NEAR(kin.getVelocity().x, 2.0, 1e-9);           // contacts see this velocity
    }

    // ---- RK4: exact-ish orbit energy vs Euler (harmonic oscillator) ----
    {
        Vector3 x(1, 0, 0), v(0, 0, 0);
        auto acc = [](const Vector3& px, const Vector3&, real) { return px * -1.0; };   // x'' = -x
        real t = 0, dt = 0.05;
        for (int i = 0; i < 2513; i++) { integrateRK4(x, v, t, dt, acc); t += dt; }     // ~20 periods
        real E = 0.5 * (v.squareMagnitude() + x.squareMagnitude());
        CHECK_NEAR(E, 0.5, 1e-3);                             // energy conserved to 0.1%
    }

    // ---- raycast vehicle drives forward and steers ----
    {
        RigidBody chassis; chassis.setMass(1200);
        Matrix3 I; I.setBlockInertiaTensor(Vector3(1.0, 0.4, 2.2), 1200); chassis.setInertiaTensor(I);
        chassis.setPosition(0, 0.75, 0); chassis.setAcceleration(0, -9.81, 0);
        chassis.setDamping(0.995, 0.8); chassis.setCanSleep(false); chassis.calculateDerivedData();
        RaycastVehicle car; car.chassis = &chassis;
        car.addWheel(Vector3(-0.85, -0.1, 1.4), true, false);
        car.addWheel(Vector3(0.85, -0.1, 1.4), true, false);
        car.addWheel(Vector3(-0.85, -0.1, -1.4), false, true);
        car.addWheel(Vector3(0.85, -0.1, -1.4), false, true);
        car.engineForce = 3000;
        for (int i = 0; i < 600; i++) { car.update(1.0 / 120); chassis.integrate(1.0 / 120); }
        CHECK(chassis.getPosition().z > 3.0);                 // accelerated forward
        CHECK(chassis.getPosition().y > 0.4 && chassis.getPosition().y < 1.1);   // suspension holds it up
        real zStraight = chassis.getPosition().z;
        car.steer = 0.35;                                     // now turn
        for (int i = 0; i < 600; i++) { car.update(1.0 / 120); chassis.integrate(1.0 / 120); }
        CHECK(real_abs(chassis.getPosition().x) > 0.5);       // path curved
        (void)zStraight;
    }

    // ---- character controller: slides along a wall, lands on ground ----
    {
        CharacterController cc;
        cc.position = Vector3(0, 1.0, 0);
        cc.planes.push_back({Vector3(0, 1, 0), 0});           // ground
        cc.planes.push_back({Vector3(-1, 0, 0), -2});         // wall at x=2 facing -x
        for (int i = 0; i < 240; i++) cc.move(Vector3(2, 0, 1), 1.0 / 60);   // run diagonally into the wall
        CHECK(cc.grounded);
        CHECK(cc.position.x < 2.01 - cc.radius + 0.02);       // blocked by the wall
        CHECK(cc.position.z > 2.0);                            // but slid along it
        CHECK_NEAR(cc.position.y, cc.halfHeight + cc.radius, 0.05);   // standing on the floor
        real y0 = cc.position.y;
        cc.jump(5);
        for (int i = 0; i < 30; i++) cc.move(Vector3(0, 0, 0), 1.0 / 60);
        CHECK(cc.position.y > y0 + 0.5);                       // airborne
        for (int i = 0; i < 200; i++) cc.move(Vector3(0, 0, 0), 1.0 / 60);
        CHECK(cc.grounded);                                    // landed again
    }

    return test::report("features");
}
