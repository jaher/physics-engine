#include "phys/phys.h"
#include "check.h"
using namespace phys;

// Contact generators that collide a primitive with the ground half-space.
struct BoxPlaneGen : ContactGenerator {
    CollisionBox* box; CollisionPlane plane; real friction, restitution;
    unsigned addContact(Contact* c, unsigned limit) const override {
        box->calculateInternals();
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = friction; d.restitution = restitution;
        CollisionDetector::boxAndHalfSpace(*box, plane, &d);
        return d.contactCount;
    }
};
struct SpherePlaneGen : ContactGenerator {
    CollisionSphere* sphere; CollisionPlane plane; real friction, restitution;
    unsigned addContact(Contact* c, unsigned limit) const override {
        sphere->calculateInternals();
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = friction; d.restitution = restitution;
        CollisionDetector::sphereAndHalfSpace(*sphere, plane, &d);
        return d.contactCount;
    }
};

static void setupBox(RigidBody& b, Vector3 pos, real mass = 1, Vector3 half = Vector3(0.5, 0.5, 0.5)) {
    b.setMass(mass); Matrix3 I; I.setBlockInertiaTensor(half, mass); b.setInertiaTensor(I);
    b.setPosition(pos); b.setAcceleration(0, -9.81, 0); b.calculateDerivedData();
}

int main() {
    CollisionPlane ground; ground.direction = Vector3(0, 1, 0); ground.offset = 0;

    // A) A box dropped onto the ground settles at rest at the right height.
    {
        RigidBody b; setupBox(b, Vector3(0, 2, 0));
        CollisionBox box; box.body = &b; box.halfSize = Vector3(0.5, 0.5, 0.5);
        BoxPlaneGen gen; gen.box = &box; gen.plane = ground; gen.friction = 0.5; gen.restitution = 0.0;
        World world(20);
        world.getRigidBodies().push_back(&b);
        world.getContactGenerators().push_back(&gen);
        for (int i = 0; i < 3000; i++) { world.startFrame(); world.runPhysics(0.005); }
        CHECK_NEAR(b.getPosition().y, 0.5, 0.05);
        CHECK(b.getVelocity().magnitude() < 0.1);
        CHECK(b.getPosition().y > 0.4);                     // never sank through the floor
    }

    // B) A bouncy box rebounds to a lower (positive) height.
    {
        RigidBody b; setupBox(b, Vector3(0, 3, 0)); b.setDamping(1, 1); b.setCanSleep(false);
        CollisionBox box; box.body = &b; box.halfSize = Vector3(0.5, 0.5, 0.5);
        BoxPlaneGen gen; gen.box = &box; gen.plane = ground; gen.friction = 0.2; gen.restitution = 0.8;
        World world(20);
        world.getRigidBodies().push_back(&b);
        world.getContactGenerators().push_back(&gen);
        bool bounced = false; real peak = 0;
        for (int i = 0; i < 3000; i++) { world.startFrame(); world.runPhysics(0.004);
            if (!bounced && b.getPosition().y < 0.6 && b.getVelocity().y > 0.5) bounced = true;
            if (bounced) peak = std::max(peak, b.getPosition().y); }
        CHECK(bounced);
        CHECK(peak > 1.2 && peak < 2.8);
    }

    // C) A sphere settles at rest on the ground.
    {
        RigidBody b; b.setMass(1); Matrix3 I; I.setInertiaTensorCoeffs(0.4, 0.4, 0.4); b.setInertiaTensor(I);
        b.setPosition(0, 2, 0); b.setAcceleration(0, -9.81, 0); b.calculateDerivedData();
        CollisionSphere sph; sph.body = &b; sph.radius = 0.5;
        SpherePlaneGen gen; gen.sphere = &sph; gen.plane = ground; gen.friction = 0.5; gen.restitution = 0.0;
        World world(10);
        world.getRigidBodies().push_back(&b);
        world.getContactGenerators().push_back(&gen);
        for (int i = 0; i < 3000; i++) { world.startFrame(); world.runPhysics(0.005); }
        CHECK_NEAR(b.getPosition().y, 0.5, 0.05);
    }

    // D) Friction slows a box sliding along the ground.
    {
        RigidBody b; setupBox(b, Vector3(0, 0.5, 0)); b.setCanSleep(false);
        b.setVelocity(4, 0, 0);
        CollisionBox box; box.body = &b; box.halfSize = Vector3(0.5, 0.5, 0.5);
        BoxPlaneGen gen; gen.box = &box; gen.plane = ground; gen.friction = 0.9; gen.restitution = 0.0;
        World world(20);
        world.getRigidBodies().push_back(&b);
        world.getContactGenerators().push_back(&gen);
        real startSpeed = 4.0;
        for (int i = 0; i < 1000; i++) { world.startFrame(); world.runPhysics(0.005); }
        CHECK(std::fabs(b.getVelocity().x) < startSpeed - 0.5);   // measurably decelerated
        CHECK(b.getPosition().y > 0.4);
    }

    // E) An almost-elastic bounce roughly conserves height (energy).
    {
        RigidBody b; setupBox(b, Vector3(0, 2, 0)); b.setDamping(1, 1); b.setCanSleep(false);
        CollisionBox box; box.body = &b; box.halfSize = Vector3(0.5, 0.5, 0.5);
        BoxPlaneGen gen; gen.box = &box; gen.plane = ground; gen.friction = 0.0; gen.restitution = 1.0;
        World world(20);
        world.getRigidBodies().push_back(&b);
        world.getContactGenerators().push_back(&gen);
        bool bounced = false; real peak = 0;
        for (int i = 0; i < 3000; i++) { world.startFrame(); world.runPhysics(0.003);
            if (!bounced && b.getPosition().y < 0.6 && b.getVelocity().y > 0.5) bounced = true;
            if (bounced) peak = std::max(peak, b.getPosition().y); }
        CHECK(peak > 1.7);                                 // near the 2.0 drop height
    }

    return test::report("resolution");
}
