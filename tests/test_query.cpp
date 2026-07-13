// Scene-query, broad-phase and serialization tests:
//   A) a fast swept sphere catches a box it would tunnel past in one discrete step;
//      a swept sphere that misses reports no hit.
//   B) Sweep-and-Prune and the dynamic AABB tree each return exactly the same
//      overlapping-pair set as an O(n^2) brute-force check (including after removals).
//   C) snapshot -> mutate -> restore returns bodies to the exact saved state.
#include "phys/query.h"
#include "phys/broadphase2.h"
#include "phys/serialize.h"
#include "check.h"
#include <random>
#include <set>
#include <vector>
using namespace phys;

static RigidBody makeBodyAt(const Vector3& p) {
    RigidBody b; b.setPosition(p); b.setOrientation(1, 0, 0, 0); b.calculateDerivedData(); return b;
}

int main() {
    // ---- A) swept-sphere continuous query ------------------------------------
    {
        RigidBody boxBody = makeBodyAt(Vector3(0, 0, 0));
        CollisionBox box; box.body = &boxBody; box.halfSize = Vector3(1, 1, 1); box.calculateInternals();
        QueryWorld qw; qw.boxes.push_back(&box);

        // A fast sphere crossing straight through the box in a single step.
        Vector3 from(-5, 0, 0), to(5, 0, 0); real r = 0.5;
        // A discrete test at both endpoints sees nothing — it tunnels.
        CHECK(!aabbOfBox(box).overlapsSphere(from, r));
        CHECK(!aabbOfBox(box).overlapsSphere(to,   r));
        // The sweep catches it with a time-of-impact strictly inside (0,1).
        SweepHit h = qw.sweepSphere(from, to, r);
        CHECK(h.hit);
        CHECK(h.toi > 0 && h.toi < 1);
        CHECK_NEAR(h.toi, 0.35, 1e-3);          // surface at x=-1, centre stops at x=-1.5
        CHECK(h.body == &boxBody);
        CHECK_NEAR(h.normal.x, -1.0, 1e-6);     // hit the -x face

        // A parallel sweep that passes above the box: no hit.
        SweepHit miss = qw.sweepSphere(Vector3(-5, 3, 0), Vector3(5, 3, 0), r);
        CHECK(!miss.hit);

        // Swept sphere vs static sphere (analytic path) and vs plane.
        RigidBody sphBody = makeBodyAt(Vector3(0, 0, 0));
        CollisionSphere s; s.body = &sphBody; s.radius = 1.0; s.calculateInternals();
        QueryWorld qs; qs.spheres.push_back(&s);
        SweepHit hs = qs.sweepSphere(Vector3(-10, 0, 0), Vector3(10, 0, 0), 0.5);
        CHECK(hs.hit); CHECK(hs.toi > 0 && hs.toi < 1);
        CHECK_NEAR(hs.toi, 0.425, 1e-6);        // centres 1.5 apart: (10-1.5)/20
        SweepHit sm = qs.sweepSphere(Vector3(-10, 5, 0), Vector3(10, 5, 0), 0.5);
        CHECK(!sm.hit);

        // Swept AABB vs static AABB: a fast box tunnelling past another.
        real t; Vector3 n;
        Aabb mover = Aabb::fromCentreHalf(Vector3(-5, 0, 0), Vector3(0.5, 0.5, 0.5));
        Aabb stat  = Aabb::fromCentreHalf(Vector3(0, 0, 0),  Vector3(1, 1, 1));
        CHECK(sweptAabbVsAabb(mover, Vector3(10, 0, 0), stat, t, n));
        CHECK(t > 0 && t < 1);
        CHECK_NEAR(t, 0.35, 1e-9);              // faces meet at -1.5 → (−1.5+5)/10
        Aabb moverMiss = Aabb::fromCentreHalf(Vector3(-5, 3, 0), Vector3(0.5, 0.5, 0.5));
        CHECK(!sweptAabbVsAabb(moverMiss, Vector3(10, 0, 0), stat, t, n));

        // Speculative contact: a fast body still a gap away from a plane, predicted
        // to hit this step, yields a contact before it penetrates.
        RigidBody fast = makeBodyAt(Vector3(0, 2, 0));
        fast.setVelocity(Vector3(0, -100, 0));
        SpeculativeContact sc = speculativeSpherePlane(&fast, 0.5,
                                                       Vector3(0, 1, 0), 0.0, 1.0 / 60.0);
        CHECK(sc.separation > 0);               // not yet touching (still a 1.5 gap)
        CHECK(sc.willHit);                      // but 100 m/s closes it within the step
        CHECK(sc.toi > 0 && sc.toi < 1);
        Contact c = sc.toContact();
        CHECK(c.penetration == 0);              // speculative: no real overlap yet
        CHECK(c.body[0] == &fast);
        // A slow body far from the plane is not predicted to hit.
        RigidBody slow = makeBodyAt(Vector3(0, 2, 0)); slow.setVelocity(Vector3(0, -1, 0));
        CHECK(!speculativeSpherePlane(&slow, 0.5, Vector3(0, 1, 0), 0.0, 1.0 / 60.0).willHit);
    }

    // ---- overlap queries ------------------------------------------------------
    {
        RigidBody b0 = makeBodyAt(Vector3(0, 0, 0));
        RigidBody b1 = makeBodyAt(Vector3(3, 0, 0));
        RigidBody b2 = makeBodyAt(Vector3(20, 0, 0));
        CollisionSphere s0; s0.body = &b0; s0.radius = 1; s0.calculateInternals();
        CollisionBox     bx; bx.body = &b1; bx.halfSize = Vector3(1, 1, 1); bx.calculateInternals();
        CollisionSphere s2; s2.body = &b2; s2.radius = 1; s2.calculateInternals();
        QueryWorld qw; qw.spheres = {&s0, &s2}; qw.boxes = {&bx};

        auto near = qw.overlapSphere(Vector3(1.5, 0, 0), 1.0);   // reaches s0 (dist 1.5≤2) & bx (dist to box .5≤1)
        CHECK(near.size() == 2);
        auto all = qw.overlapAabb(Aabb::fromCentreHalf(Vector3(2, 0, 0), Vector3(10, 10, 10)));
        CHECK(all.size() == 2);                                  // s0 and bx, not the far s2
    }

    // ---- B) broad phase: SAP & AABB-tree vs brute force -----------------------
    {
        std::mt19937 rng(12345);
        std::uniform_real_distribution<double> pos(-10, 10), half(0.2, 1.5);
        const int N = 220;
        std::vector<Aabb> boxes;
        for (int i = 0; i < N; i++) {
            Vector3 c(pos(rng), pos(rng), pos(rng));
            Vector3 hs(half(rng), half(rng), half(rng));
            boxes.push_back(Aabb::fromCentreHalf(c, hs));
        }
        // brute force O(n^2)
        std::set<BroadPair> brute;
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++)
                if (boxes[i].overlaps(boxes[j])) brute.insert(makePair(i, j));

        // sweep and prune
        SweepAndPrune sap;
        for (int i = 0; i < N; i++) sap.add(boxes[i], i);
        auto sapPairs = sap.overlappingPairs();
        std::set<BroadPair> sapSet(sapPairs.begin(), sapPairs.end());
        CHECK(sapSet == brute);

        // dynamic AABB tree
        DynamicAabbTree tree(0.0);
        std::vector<int> handle(N);
        for (int i = 0; i < N; i++) handle[i] = tree.insert(i, boxes[i]);
        auto treePairs = tree.queryPairs();
        std::set<BroadPair> treeSet(treePairs.begin(), treePairs.end());
        CHECK(treeSet == brute);
        CHECK(tree.leafCount() == N);

        // queryAABB agrees with a linear scan
        Aabb probe = Aabb::fromCentreHalf(Vector3(1, -2, 3), Vector3(4, 4, 4));
        std::set<int> treeHits, linHits;
        for (int id : tree.queryAABB(probe)) treeHits.insert(id);
        for (int i = 0; i < N; i++) if (boxes[i].overlaps(probe)) linHits.insert(i);
        CHECK(treeHits == linHits);

        // remove ~half, then the surviving pair set must still match brute force
        std::set<int> removed;
        for (int i = 0; i < N; i += 2) { tree.remove(handle[i]); removed.insert(i); }
        std::set<BroadPair> brute2;
        for (int i = 0; i < N; i++) {
            if (removed.count(i)) continue;
            for (int j = i + 1; j < N; j++) {
                if (removed.count(j)) continue;
                if (boxes[i].overlaps(boxes[j])) brute2.insert(makePair(i, j));
            }
        }
        auto treePairs2 = tree.queryPairs();
        std::set<BroadPair> treeSet2(treePairs2.begin(), treePairs2.end());
        CHECK(treeSet2 == brute2);
        CHECK(tree.leafCount() == N - (int)removed.size());
    }

    // ---- C) serialize -> mutate -> restore ------------------------------------
    {
        RigidBody b0, b1;
        b0.setMass(2); b1.setMass(3);
        Matrix3 it; it.setBlockInertiaTensor(Vector3(1, 1, 1), 2);
        b0.setInertiaTensor(it); b1.setInertiaTensor(it);
        b0.setPosition(1.5, -2.25, 3.125);   b1.setPosition(-7.0, 0.5, 11.0);
        b0.setVelocity(0.5, -0.25, 0.125);   b1.setVelocity(-3.0, 4.0, -0.5);
        b0.setRotation(0.1, 0.2, -0.3);      b1.setRotation(-0.4, 0.05, 0.6);
        b0.setOrientation(Quaternion(0.5, 0.5, 0.5, 0.5));           // already unit
        b1.setOrientation(Quaternion(0.9238795, 0.0, 0.3826834, 0)); // ~45° about y
        b0.setAcceleration(Vector3::GRAVITY); b1.setAcceleration(Vector3::GRAVITY);
        b0.calculateDerivedData(); b1.calculateDerivedData();

        std::vector<RigidBody*> bodies = {&b0, &b1};

        // capture the exact saved state
        struct S { Vector3 p, v, w; Quaternion q; } saved[2];
        for (int i = 0; i < 2; i++) {
            saved[i] = {bodies[i]->getPosition(), bodies[i]->getVelocity(),
                        bodies[i]->getRotation(), bodies[i]->getOrientation()};
        }

        std::vector<unsigned char> buf = Snapshot::snapshot(bodies);
        CHECK(buf.size() == Snapshot::byteSize(2));

        // mutate: run the sim forward, scrambling every field
        for (int s = 0; s < 240; s++) { b0.integrate(1.0 / 60); b1.integrate(1.0 / 60); }
        CHECK(bodies[0]->getPosition().y != saved[0].p.y);   // it really moved

        // restore
        CHECK(Snapshot::restore(buf, bodies));
        for (int i = 0; i < 2; i++) {
            Vector3 p = bodies[i]->getPosition(), v = bodies[i]->getVelocity(), w = bodies[i]->getRotation();
            Quaternion q = bodies[i]->getOrientation();
            // position & velocities are restored bitwise-exactly
            CHECK(p.x == saved[i].p.x && p.y == saved[i].p.y && p.z == saved[i].p.z);
            CHECK(v.x == saved[i].v.x && v.y == saved[i].v.y && v.z == saved[i].v.z);
            CHECK(w.x == saved[i].w.x && w.y == saved[i].w.y && w.z == saved[i].w.z);
            // orientation is a unit quaternion reproduced to full double precision
            CHECK_NEAR(q.r, saved[i].q.r, 1e-12);
            CHECK_NEAR(q.i, saved[i].q.i, 1e-12);
            CHECK_NEAR(q.j, saved[i].q.j, 1e-12);
            CHECK_NEAR(q.k, saved[i].q.k, 1e-12);
        }

        // a size mismatch or bad magic is rejected
        std::vector<RigidBody*> one = {&b0};
        CHECK(!Snapshot::restore(buf, one));
        std::vector<unsigned char> corrupt = buf; corrupt[0] ^= 0xFF;
        CHECK(!Snapshot::restore(corrupt, bodies));
    }

    return test::report("query");
}
