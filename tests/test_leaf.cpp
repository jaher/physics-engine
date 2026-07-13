// Falling burning leaf: the oak silhouette is a proper shape, the tumbling-plate
// aerodynamics give a bounded fluttering descent, and combustion consumes the leaf.
#include "phys/leaf.h"
#include "check.h"
#include <cmath>
using namespace phys;

int main() {
    // A) oak mask is a shape (not the whole rectangle), veins exist
    {
        FallingLeaf L; L.build(0.12, 0.08, 42);
        CHECK(L.maskCount > 0);
        CHECK(L.maskCount < L.GW * L.GH);                 // a silhouette, not a filled rect
        CHECK(L.leafArea > 0);
        real vsum = 0; for (float v : L.vein) vsum += v; CHECK(vsum > 0);
    }

    // B) aerodynamics: the leaf falls at a bounded speed, flutters sideways, and its
    //    orientation frame stays orthonormal (no blow-up, no NaN)
    {
        FallingLeaf L; L.build(0.12, 0.08, 7);
        L.pos = Vector3(0, 5, 0); L.vel = Vector3(0.1, 0, 0);
        L.e2 = Vector3(0.2, 0.95, 0.2); L.e2.normalise(); L.e0 = (L.e2 % Vector3(0, 1, 0.3)).unit(); L.e1 = (L.e2 % L.e0);
        real maxSpeed = 0; real horiz = 0;
        for (int s = 0; s < 1440; s++) { L.step(1.0 / 240, Vector3()); maxSpeed = std::max(maxSpeed, (real)std::sqrt(L.vel.squareMagnitude())); }
        horiz = std::sqrt(L.pos.x * L.pos.x + L.pos.z * L.pos.z);
        CHECK(std::isfinite(L.pos.y) && std::isfinite(maxSpeed));
        CHECK(L.pos.y < 5.0);                             // it fell
        CHECK(maxSpeed < 4.0);                            // drag keeps it slow (not free-fall ~19 m/s over 6 s)
        CHECK(horiz > 0.05);                              // it drifted sideways (fluttered), not straight down
        CHECK_NEAR(L.e0.magnitude(), 1.0, 1e-5);          // frame stays orthonormal
        CHECK(std::fabs(L.e0 * L.e1) < 1e-5);
        CHECK_NEAR((L.e0 % L.e1).magnitude(), 1.0, 1e-5);
    }

    // C) an un-ignited leaf never loses fuel
    {
        FallingLeaf L; L.build(0.12, 0.08, 3); L.pos = Vector3(0, 5, 0);
        for (int s = 0; s < 480; s++) L.step(1.0 / 240, Vector3());
        CHECK(L.aliveFrac() == 1.0);
    }

    // D) an ignited leaf burns: fuel is consumed and it is eventually consumed away
    {
        FallingLeaf L; L.build(0.12, 0.08, 9); L.pos = Vector3(0, 5, 0);
        L.burnRate = 0.5; L.spread = 1.6; L.igniteAt = 0.0;
        real start = L.aliveFrac();
        for (int s = 0; s < 240; s++) L.step(1.0 / 240, Vector3());   // 1 s in
        CHECK(L.ignited);
        CHECK(L.aliveFrac() < start);                                 // burning eats the leaf
        for (int s = 0; s < 2400; s++) L.step(1.0 / 240, Vector3());  // let it finish
        CHECK(L.consumed());
    }

    return test::report("leaf");
}
