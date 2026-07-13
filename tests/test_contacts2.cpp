// Persistent manifold + warm-starting, material combine modes, and rolling /
// spinning friction — the Bullet/PhysX-quality contact-solver upgrade.
#include "check.h"
#include "phys/body.h"
#include "phys/material.h"
#include "phys/contacts.h"
#include "phys/collide_fine.h"
#include <cmath>
#include <cstdio>
using namespace phys;

// ---- a resting stack of three unit boxes on the ground half-space -----------
// Driven directly against a ContactResolver (not World) so we control the exact
// iteration budget and can turn warm-starting on/off. Runs the stack to rest with a
// generous, stable budget, then reports (out) the total residual penetration
// (box0↓ground + box0/box1 overlap + box1/box2 overlap) and the number of velocity
// solver iterations the final, settled frame needed — the warm-start payoff.
static void runStack(bool warm, int frames, real& penOut, unsigned& lastVelItersOut) {
    CollisionPlane ground; ground.direction = Vector3(0, 1, 0); ground.offset = 0;
    RigidBody b[3]; CollisionBox cb[3];
    for (int i = 0; i < 3; i++) {
        b[i].setMass(1);
        Matrix3 I; I.setBlockInertiaTensor(Vector3(0.5, 0.5, 0.5), 1); b[i].setInertiaTensor(I);
        b[i].setPosition(0, 0.5 + i * 1.0, 0);
        b[i].setAcceleration(0, -9.81, 0);
        b[i].setCanSleep(false);
        b[i].calculateDerivedData();
        cb[i].body = &b[i]; cb[i].halfSize = Vector3(0.5, 0.5, 0.5);
    }
    ContactResolver resolver(64u, 64u);        // generous budget: both configs stay stable
    resolver.setWarmStarting(warm);
    Contact contacts[32];
    const real dt = 0.01;
    for (int f = 0; f < frames; f++) {
        for (int i = 0; i < 3; i++) b[i].integrate(dt);
        for (int i = 0; i < 3; i++) cb[i].calculateInternals();
        CollisionData d; d.contactArray = contacts; d.reset(32); d.friction = 0.5; d.restitution = 0;
        CollisionDetector::boxAndHalfSpace(cb[0], ground, &d);
        CollisionDetector::boxAndBox(cb[0], cb[1], &d);
        CollisionDetector::boxAndBox(cb[1], cb[2], &d);
        unsigned nc = d.contactCount;
        if (nc) resolver.resolveContacts(contacts, nc, dt);
    }
    penOut = 0;
    penOut += std::max((real)0, (real)0.5 - b[0].getPosition().y);
    penOut += std::max((real)0, (real)1.0 - (b[1].getPosition().y - b[0].getPosition().y));
    penOut += std::max((real)0, (real)1.0 - (b[2].getPosition().y - b[1].getPosition().y));
    lastVelItersOut = resolver.velocityIterationsUsed;   // iterations spent on the settled frame
}

// ---- a pure-rolling sphere on the ground; returns its final |angular velocity| --
static real rollingSphereFinalOmega(real rollFric) {
    CollisionPlane ground; ground.direction = Vector3(0, 1, 0); ground.offset = 0;
    RigidBody b; b.setMass(1);
    Matrix3 I; I.setInertiaTensorCoeffs(0.4, 0.4, 0.4); b.setInertiaTensor(I);
    b.setPosition(0, 0.5, 0);
    b.setAcceleration(0, -9.81, 0);
    b.setDamping(1, 1);            // isolate rolling friction: no velocity damping
    b.setCanSleep(false);          // must not "stop" by falling asleep
    const real w0 = 4.0, r = 0.5;
    b.setRotation(0, 0, -w0);      // spin about z ...
    b.setVelocity(w0 * r, 0, 0);   // ... matched to +x translation → pure rolling (no slip)
    b.calculateDerivedData();
    CollisionSphere s; s.body = &b; s.radius = r;
    ContactResolver resolver(10u, 10u);
    Contact contacts[4];
    const real dt = 0.01;
    for (int f = 0; f < 4000; f++) {
        b.integrate(dt);
        s.calculateInternals();
        CollisionData d; d.contactArray = contacts; d.reset(4); d.friction = 0.6; d.restitution = 0;
        CollisionDetector::sphereAndHalfSpace(s, ground, &d);
        unsigned nc = d.contactCount;
        if (nc) { contacts[0].setRollingFriction(rollFric, 0); resolver.resolveContacts(contacts, nc, dt); }
    }
    return b.getRotation().magnitude();
}

int main() {
    // (1) Material combine modes return the expected numbers.
    {
        Material a(0.4, 0.2, 0.10, 0.05);
        Material b(0.6, 0.8, 0.30, 0.15);

        Material avg = combine(a, b);                       // default = Average
        CHECK_NEAR(avg.friction,        0.5,  1e-12);
        CHECK_NEAR(avg.restitution,     0.5,  1e-12);
        CHECK_NEAR(avg.rollingFriction, 0.2,  1e-12);
        CHECK_NEAR(avg.spinFriction,    0.10, 1e-12);

        Material mn = combine(a, b, CombineMode::Min);
        CHECK_NEAR(mn.friction,    0.4, 1e-12);
        CHECK_NEAR(mn.restitution, 0.2, 1e-12);

        Material mx = combine(a, b, CombineMode::Max);
        CHECK_NEAR(mx.friction,    0.6, 1e-12);
        CHECK_NEAR(mx.restitution, 0.8, 1e-12);

        Material ml = combine(a, b, CombineMode::Multiply);
        CHECK_NEAR(ml.friction,        0.24, 1e-12);
        CHECK_NEAR(ml.restitution,     0.16, 1e-12);
        CHECK_NEAR(ml.rollingFriction, 0.03, 1e-12);
    }

    // (2) Warm-starting: a settled stack keeps ≤ the penetration of a cold solve
    //     and reaches rest in strictly fewer velocity-solver iterations per frame.
    {
        real coldPen = 0, warmPen = 0; unsigned coldIt = 0, warmIt = 0;
        runStack(false, 400, coldPen, coldIt);
        runStack(true,  400, warmPen, warmIt);
        std::printf("  [warm-start] settled stack: cold pen=%.6f iters=%u   warm pen=%.6f iters=%u\n",
                    coldPen, coldIt, warmPen, warmIt);
        CHECK(warmPen <= coldPen + 1e-9);   // never worse penetration than the cold solve
        CHECK(coldIt > 0);                  // the cold solve really does re-iterate every frame
        CHECK(warmIt < coldIt);             // warm-started stack settles in fewer iterations
    }

    // (3) Rolling friction brings a rolling sphere to rest; zero rolling friction
    //     leaves it rolling.
    {
        real rolling = rollingSphereFinalOmega(0.0);
        real stopped = rollingSphereFinalOmega(0.3);
        std::printf("  [rolling]   final |omega|: no-friction=%.4f  rolling-friction=%.4f\n", rolling, stopped);
        CHECK(rolling > 2.0);               // still rolling without rolling friction
        CHECK(stopped < 0.2);               // slowed to rest with rolling friction
        CHECK(stopped < rolling);
    }

    return test::report("contacts2");
}
