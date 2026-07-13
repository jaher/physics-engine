#include "phys/phys.h"
#include "check.h"
using namespace phys;

struct BoxBoxGen : ContactGenerator {
    CollisionBox* a; CollisionBox* b; real friction, restitution;
    unsigned addContact(Contact* c, unsigned limit) const override {
        a->calculateInternals(); b->calculateInternals();
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = friction; d.restitution = restitution;
        CollisionDetector::boxAndBox(*a, *b, &d);
        return d.contactCount;
    }
};
static void setupBox(RigidBody& b, Vector3 pos, Vector3 accel = Vector3(0, 0, 0)) {
    b.setMass(1); Matrix3 I; I.setBlockInertiaTensor(Vector3(0.5, 0.5, 0.5), 1); b.setInertiaTensor(I);
    b.setPosition(pos); b.setAcceleration(accel); b.setCanSleep(false); b.calculateDerivedData();
}

int main() {
    // A) Two interpenetrating boxes are pushed apart to (roughly) touching.
    {
        RigidBody A, B; setupBox(A, Vector3(0, 0, 0)); setupBox(B, Vector3(0.6, 0, 0));
        CollisionBox ba; ba.body = &A; ba.halfSize = Vector3(0.5, 0.5, 0.5);
        CollisionBox bb; bb.body = &B; bb.halfSize = Vector3(0.5, 0.5, 0.5);
        BoxBoxGen gen; gen.a = &ba; gen.b = &bb; gen.friction = 0; gen.restitution = 0;
        World world(20);
        world.getRigidBodies().push_back(&A); world.getRigidBodies().push_back(&B);
        world.getContactGenerators().push_back(&gen);
        for (int i = 0; i < 300; i++) { world.startFrame(); world.runPhysics(0.01); }
        real dx = std::fabs(A.getPosition().x - B.getPosition().x);
        CHECK(dx > 0.9 && dx < 1.15);                       // separated to ~1.0
    }

    // B) A box dropped onto a fixed box comes to rest on top of it (box-box).
    {
        RigidBody bottom; setupBox(bottom, Vector3(0, 0.5, 0));
        bottom.setInverseMass(0); bottom.setInverseInertiaTensor(Matrix3()); bottom.calculateDerivedData();
        RigidBody top; setupBox(top, Vector3(0, 3, 0), Vector3(0, -9.81, 0));
        CollisionBox bb; bb.body = &bottom; bb.halfSize = Vector3(0.5, 0.5, 0.5);
        CollisionBox bt; bt.body = &top; bt.halfSize = Vector3(0.5, 0.5, 0.5);
        BoxBoxGen gen; gen.a = &bb; gen.b = &bt; gen.friction = 0.4; gen.restitution = 0;
        World world(20);
        world.getRigidBodies().push_back(&bottom); world.getRigidBodies().push_back(&top);
        world.getContactGenerators().push_back(&gen);
        for (int i = 0; i < 4000; i++) { world.startFrame(); world.runPhysics(0.004); }
        CHECK(top.getPosition().y > 1.2 && top.getPosition().y < 1.7);   // resting on the lower box
        CHECK(bottom.getPosition().y == 0.5);                            // fixed body never moved
    }

    return test::report("stacking");
}
