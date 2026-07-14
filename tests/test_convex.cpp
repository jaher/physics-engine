// Convex / mesh geometry: GJK boolean + distance, EPA penetration, cylinder &
// cone vs plane, and a sphere resting on a static triangle-mesh ramp.
#include "phys/collide_convex.h"
#include "check.h"
#include <cmath>
using namespace phys;

static bool vnear(const Vector3& a, const Vector3& b, double t = 1e-6) {
    return test::near(a.x, b.x, t) && test::near(a.y, b.y, t) && test::near(a.z, b.z, t);
}
static RigidBody* place(RigidBody& b, Vector3 pos, Quaternion q = Quaternion()) {
    b.setMass(1); Matrix3 I; I.setBlockInertiaTensor(Vector3(1, 1, 1), 1); b.setInertiaTensor(I);
    b.setPosition(pos); b.setOrientation(q); b.calculateDerivedData(); return &b;
}

int main() {
    Contact contacts[64];
    CollisionData data; data.contactArray = contacts;

    // A) two overlapping convex hulls (unit boxes overlapping 0.5 along x): EPA
    //    reports penetration ~0.5 with a normal pointing from b toward a (−x).
    {
        RigidBody a, b; place(a, Vector3(0, 0, 0)); place(b, Vector3(1.5, 0, 0));
        CollisionConvex ca; ca.body = &a; ca.setBox(Vector3(1, 1, 1)); ca.calculateInternals();
        CollisionConvex cb; cb.body = &b; cb.setBox(Vector3(1, 1, 1)); cb.calculateInternals();
        CHECK(ConvexCollision::intersect(ca, cb));
        data.reset(64);
        unsigned n = ConvexCollision::convexAndConvex(ca, cb, &data);
        CHECK(n == 1);
        CHECK(contacts[0].penetration > 0);
        CHECK_NEAR(contacts[0].penetration, 0.5, 1e-3);
        CHECK_NEAR(std::fabs(contacts[0].contactNormal.x), 1.0, 1e-3);  // separating along x
        CHECK(contacts[0].contactNormal.x < 0);                         // from b(+x) toward a
        Vector3 sep = a.getPosition() - b.getPosition();
        CHECK(contacts[0].contactNormal * sep > 0);                     // points body[1]→body[0]
    }

    // B) two separated hulls report no contact and no intersection.
    {
        RigidBody a, b; place(a, Vector3(0, 0, 0)); place(b, Vector3(5, 0, 0));
        CollisionConvex ca; ca.body = &a; ca.setBox(Vector3(1, 1, 1)); ca.calculateInternals();
        CollisionConvex cb; cb.body = &b; cb.setBox(Vector3(1, 1, 1)); cb.calculateInternals();
        CHECK(!ConvexCollision::intersect(ca, cb));
        data.reset(64);
        unsigned n = ConvexCollision::convexAndConvex(ca, cb, &data);
        CHECK(n == 0);
    }

    // C) GJK distance for a known pair: half-0.5 boxes at x=0 and x=3 → gap 2.0.
    {
        RigidBody a, b; place(a, Vector3(0, 0, 0)); place(b, Vector3(3, 0, 0));
        CollisionConvex ca; ca.body = &a; ca.setBox(Vector3(0.5, 0.5, 0.5)); ca.calculateInternals();
        CollisionConvex cb; cb.body = &b; cb.setBox(Vector3(0.5, 0.5, 0.5)); cb.calculateInternals();
        Vector3 pA, pB;
        real d = ConvexCollision::distance(ca, cb, &pA, &pB);
        CHECK_NEAR(d, 2.0, 1e-6);
        CHECK_NEAR(pA.x, 0.5, 1e-6);       // witness on a's +x face
        CHECK_NEAR(pB.x, 2.5, 1e-6);       // witness on b's −x face
        // ...and touching hulls have zero distance.
        RigidBody c; place(c, Vector3(1.0, 0, 0));   // half-0.5 box just touching a at x=0.5
        CollisionConvex cc; cc.body = &c; cc.setBox(Vector3(0.5, 0.5, 0.5)); cc.calculateInternals();
        CHECK(ConvexCollision::distance(ca, cc) < 1e-6);
    }

    // D) sphere resting on a tilted triangle-mesh ramp: contact normal ≈ ramp normal.
    {
        const double ang = 0.3;                       // ramp tilt about the z axis
        const double ca = std::cos(ang), sa = std::sin(ang), L = 2.0;
        Vector3 nramp(-sa, ca, 0);                    // known ramp normal (plane through origin)
        // flat quad in the xz-plane rotated about z: (x,0,z) → (x·cos, x·sin, z)
        Vector3 p0(-L * ca, -L * sa, -L), p1(L * ca, L * sa, -L),
                p2( L * ca,  L * sa,  L), p3(-L * ca, -L * sa,  L);
        CollisionTriMesh mesh;
        mesh.addTriangle(p0, p1, p2);
        mesh.addTriangle(p0, p2, p3);
        mesh.build();                                 // exercise the uniform-grid midphase

        RigidBody s; place(s, nramp * 0.8);           // sphere r=1 sunk 0.2 into the ramp
        CollisionSphere sph; sph.body = &s; sph.radius = 1; sph.calculateInternals();
        data.reset(64);
        unsigned n = ConvexCollision::sphereAndTriMesh(sph, mesh, &data);
        CHECK(n == 1);
        CHECK(vnear(contacts[0].contactNormal, nramp, 1e-6));
        CHECK_NEAR(contacts[0].penetration, 0.2, 1e-6);
        CHECK(vnear(contacts[0].contactPoint, Vector3(), 1e-6));   // touches at the ramp centre

        // a sphere hovering above the ramp reports nothing.
        RigidBody s2; place(s2, nramp * 2.0);
        CollisionSphere sph2; sph2.body = &s2; sph2.radius = 1; sph2.calculateInternals();
        data.reset(64);
        CHECK(ConvexCollision::sphereAndTriMesh(sph2, mesh, &data) == 0);
    }

    // E) cylinder penetrating a ground half-space: a contact with an upward normal.
    {
        RigidBody rb; place(rb, Vector3(0, 0.9, 0));      // bottom cap sunk 0.1 below y=0
        CollisionCylinder cyl; cyl.body = &rb; cyl.radius = 0.5; cyl.halfHeight = 1;
        cyl.calculateInternals();
        CollisionPlane plane; plane.direction = Vector3(0, 1, 0); plane.offset = 0;
        data.reset(64);
        unsigned n = ConvexCollision::cylinderAndHalfSpace(cyl, plane, &data);
        CHECK(n >= 1);
        CHECK(vnear(contacts[0].contactNormal, Vector3(0, 1, 0)));
        CHECK_NEAR(contacts[0].penetration, 0.1, 1e-9);
    }

    // F) cone resting base-down, penetrating the ground: upward-normal contact.
    {
        RigidBody rb; place(rb, Vector3(0, 0.9, 0));      // base sunk 0.1 below y=0, apex up
        CollisionCone cone; cone.body = &rb; cone.radius = 0.5; cone.halfHeight = 1;
        cone.calculateInternals();
        CollisionPlane plane; plane.direction = Vector3(0, 1, 0); plane.offset = 0;
        data.reset(64);
        unsigned n = ConvexCollision::coneAndHalfSpace(cone, plane, &data);
        CHECK(n >= 1);
        CHECK(vnear(contacts[0].contactNormal, Vector3(0, 1, 0)));
        CHECK_NEAR(contacts[0].penetration, 0.1, 1e-9);
    }

    // G) a box dropped onto the flat mesh floor generates corner contacts.
    {
        CollisionTriMesh floor;                           // two triangles share the z=x diagonal
        floor.addTriangle(Vector3(-5, 0, -5), Vector3(5, 0, -5), Vector3(5, 0, 5));
        floor.addTriangle(Vector3(-5, 0, -5), Vector3(5, 0, 5), Vector3(-5, 0, 5));
        RigidBody rb; place(rb, Vector3(0, 0.4, 2));       // sits fully on the z>x triangle, off the seam
        CollisionBox box; box.body = &rb; box.halfSize = Vector3(0.5, 0.5, 0.5); box.calculateInternals();
        data.reset(64);
        unsigned n = ConvexCollision::boxAndTriMesh(box, floor, &data);
        CHECK(n == 4);                                    // four bottom corners
        for (unsigned i = 0; i < n; i++) {
            CHECK(vnear(contacts[i].contactNormal, Vector3(0, 1, 0)));
            CHECK_NEAR(contacts[i].penetration, 0.1, 1e-9);
        }
    }

    // H) cylinder vs sphere (rounded-side approximation) gives a sane contact.
    {
        RigidBody rc, rs; place(rc, Vector3(0, 0, 0)); place(rs, Vector3(0.9, 0, 0));
        CollisionCylinder cyl; cyl.body = &rc; cyl.radius = 0.5; cyl.halfHeight = 1; cyl.calculateInternals();
        CollisionSphere sph; sph.body = &rs; sph.radius = 0.5; sph.calculateInternals();
        data.reset(64);
        unsigned n = ConvexCollision::cylinderAndSphere(cyl, sph, &data);
        CHECK(n == 1);
        CHECK_NEAR(contacts[0].penetration, 0.1, 1e-6);   // 0.5 + 0.5 − 0.9
        CHECK(vnear(contacts[0].contactNormal, Vector3(-1, 0, 0), 1e-6));
    }

    // F) convex hull vs a half-space: the four bottom corners below the floor make
    //    contacts (normal up, penetration = depth); lifted clear → none.
    {
        RigidBody a; place(a, Vector3(0, 0.3, 0));
        CollisionConvex ca; ca.body = &a; ca.setBox(Vector3(1, 1, 1)); ca.calculateInternals();
        CollisionPlane floor; floor.direction = Vector3(0, 1, 0); floor.offset = 0;
        data.reset(64);
        unsigned n = ConvexCollision::convexAndHalfSpace(ca, floor, &data);
        CHECK(n == 4);                                    // 4 bottom corners at y = −0.7
        CHECK_NEAR(contacts[0].penetration, 0.7, 1e-6);
        CHECK(vnear(contacts[0].contactNormal, Vector3(0, 1, 0), 1e-6));
        place(a, Vector3(0, 2.0, 0)); ca.calculateInternals();
        data.reset(64);
        CHECK(ConvexCollision::convexAndHalfSpace(ca, floor, &data) == 0);   // lifted clear
    }

    return test::report("convex");
}
