#include "phys/phys.h"
#include "check.h"
using namespace phys;

static bool vnear(const Vector3& a, const Vector3& b, double t = 1e-6) {
    return test::near(a.x, b.x, t) && test::near(a.y, b.y, t) && test::near(a.z, b.z, t);
}
static RigidBody* place(RigidBody& b, Vector3 pos, Quaternion q = Quaternion()) {
    b.setMass(1); Matrix3 I; I.setBlockInertiaTensor(Vector3(1, 1, 1), 1); b.setInertiaTensor(I);
    b.setPosition(pos); b.setOrientation(q); b.calculateDerivedData(); return &b;
}

int main() {
    Contact contacts[16];
    CollisionData data; data.contactArray = contacts;

    // A) sphere vs half-space (ground at y=0)
    {
        RigidBody rb; place(rb, Vector3(0, 0.5, 0));
        CollisionSphere s; s.body = &rb; s.radius = 1; s.calculateInternals();
        CollisionPlane plane; plane.direction = Vector3(0, 1, 0); plane.offset = 0;
        data.reset(16); data.friction = 0; data.restitution = 0;
        unsigned n = CollisionDetector::sphereAndHalfSpace(s, plane, &data);
        CHECK(n == 1);
        CHECK_NEAR(contacts[0].penetration, 0.5, 1e-9);
        CHECK(vnear(contacts[0].contactNormal, Vector3(0, 1, 0)));
    }

    // B) sphere vs sphere
    {
        RigidBody a, b; place(a, Vector3(0, 0, 0)); place(b, Vector3(1.5, 0, 0));
        CollisionSphere sa; sa.body = &a; sa.radius = 1; sa.calculateInternals();
        CollisionSphere sb; sb.body = &b; sb.radius = 1; sb.calculateInternals();
        CHECK(IntersectionTests::sphereAndSphere(sa, sb));
        data.reset(16);
        unsigned n = CollisionDetector::sphereAndSphere(sa, sb, &data);
        CHECK(n == 1);
        CHECK_NEAR(contacts[0].penetration, 0.5, 1e-9);
        CHECK(vnear(contacts[0].contactNormal, Vector3(-1, 0, 0)));      // from b toward a
        CHECK(vnear(contacts[0].contactPoint, Vector3(0.75, 0, 0)));     // overlap midpoint
    }

    // C) box vs half-space: the four bottom corners of a unit box sunk 0.25 into ground
    {
        RigidBody rb; place(rb, Vector3(0, 0.25, 0));
        CollisionBox box; box.body = &rb; box.halfSize = Vector3(0.5, 0.5, 0.5); box.calculateInternals();
        CollisionPlane plane; plane.direction = Vector3(0, 1, 0); plane.offset = 0;
        data.reset(16);
        unsigned n = CollisionDetector::boxAndHalfSpace(box, plane, &data);
        CHECK(n == 4);
        for (unsigned i = 0; i < n; i++) { CHECK_NEAR(contacts[i].penetration, 0.25, 1e-9);
            CHECK(vnear(contacts[i].contactNormal, Vector3(0, 1, 0))); }
    }

    // D) box vs sphere touching a face
    {
        RigidBody rb, rs; place(rb, Vector3(0, 0, 0)); place(rs, Vector3(1.2, 0, 0));
        CollisionBox box; box.body = &rb; box.halfSize = Vector3(1, 1, 1); box.calculateInternals();
        CollisionSphere sph; sph.body = &rs; sph.radius = 0.5; sph.calculateInternals();
        data.reset(16);
        unsigned n = CollisionDetector::boxAndSphere(box, sph, &data);
        CHECK(n == 1);
        CHECK_NEAR(contacts[0].penetration, 0.3, 1e-6);                  // 1 + 0.5 - 1.2
        CHECK(vnear(contacts[0].contactPoint, Vector3(1, 0, 0), 1e-6));
    }

    // E) box vs box, face overlap along x
    {
        RigidBody a, b; place(a, Vector3(0, 0, 0)); place(b, Vector3(1.5, 0, 0));
        CollisionBox ba; ba.body = &a; ba.halfSize = Vector3(1, 1, 1); ba.calculateInternals();
        CollisionBox bb; bb.body = &b; bb.halfSize = Vector3(1, 1, 1); bb.calculateInternals();
        data.reset(16);
        unsigned n = CollisionDetector::boxAndBox(ba, bb, &data);
        CHECK(n == 1);
        CHECK_NEAR(contacts[0].penetration, 0.5, 1e-6);
        CHECK_NEAR(std::fabs(contacts[0].contactNormal.x), 1.0, 1e-6);
    }

    // F) BVH broad phase finds the overlapping pair (and not the far one)
    {
        RigidBody a, b, c; place(a, Vector3(0, 0, 0)); place(b, Vector3(1, 0, 0)); place(c, Vector3(50, 0, 0));
        BVHNode<BoundingSphere>* root = new BVHNode<BoundingSphere>(nullptr, BoundingSphere(a.getPosition(), 1), &a);
        root->insert(&b, BoundingSphere(b.getPosition(), 1));
        root->insert(&c, BoundingSphere(c.getPosition(), 1));
        PotentialContact pc[16];
        unsigned n = root->getPotentialContacts(pc, 16);
        CHECK(n == 1);
        CHECK((pc[0].body[0] == &a && pc[0].body[1] == &b) || (pc[0].body[0] == &b && pc[0].body[1] == &a));
        delete root;
    }

    return test::report("collision");
}
