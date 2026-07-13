// Advanced / paradigm features: forward-mode autodiff (differentiable physics),
// an island-parallel solver that is bit-identical to serial, and approximate
// convex decomposition of a concave polygon.
#include "phys/autodiff.h"
#include "phys/parallel.h"
#include "phys/decompose.h"
#include "check.h"
#include <cmath>
#include <thread>
#include <vector>
using namespace phys;

// convexity via the test's own cross-product-sign rule (independent of the header)
static bool crossSignsAgree(const Polygon2& p) {
    int n = (int)p.size(); if (n < 3) return false;
    int sign = 0;
    for (int i = 0; i < n; i++) {
        Vector3 a = p[(i - 1 + n) % n], b = p[i], c = p[(i + 1) % n];
        real cr = (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
        if (cr > 1e-9)      { if (sign < 0) return false; sign = 1; }
        else if (cr < -1e-9){ if (sign > 0) return false; sign = -1; }
    }
    return true;
}

int main() {
    // ============================ A) DIFFERENTIABLE PHYSICS =====================
    // A dual number's `.d` is the exact derivative; check it against a central
    // finite-difference gradient of the very same simulation.

    // A0) the primitive rules are right: d/dx sin = cos, d/dx (x²) via pow, etc.
    {
        real x0 = 0.7;
        Dual s = sin(variable(x0));
        CHECK_NEAR(s.v, std::sin(x0), 1e-12);
        CHECK_NEAR(s.d, std::cos(x0), 1e-12);          // d(sin)/dx = cos
        Dual q = pow(variable(x0), 3.0);
        CHECK_NEAR(q.d, 3 * x0 * x0, 1e-12);           // d(x³)/dx = 3x²
        Dual g = variable(x0) * variable(x0);          // product rule → 2x
        CHECK_NEAR(g.d, 2 * x0, 1e-12);
    }

    // A1) d(range)/d(angle) of a projectile — Dual vs central FD vs closed form.
    {
        real v0 = 20.0, g = 9.81, ang = 0.6;
        Dual r = projectileRange(variable(ang), v0, g);
        real gradAD = r.d;
        real h = 1e-6;
        real fp = projectileRange<real>(ang + h, v0, g);
        real fm = projectileRange<real>(ang - h, v0, g);
        real gradFD = (fp - fm) / (2 * h);
        real gradExact = v0 * v0 * 2 * std::cos(2 * ang) / g;   // d/dθ [v₀²sin2θ/g]
        CHECK_NEAR(gradAD, gradFD, 1e-4);              // autodiff matches finite diff
        CHECK_NEAR(gradAD, gradExact, 1e-9);           // …and the analytic gradient
        CHECK_NEAR(r.v, v0 * v0 * std::sin(2 * ang) / g, 1e-9);  // the value is right too
    }

    // A2) d(settle)/d(stiffness) of a damped spring, integrated over 400 steps —
    //     the gradient flows through every time step of the discrete simulation.
    {
        real m = 1.0, c = 0.3, x0 = 1.0, dt = 0.01; int steps = 400; real k = 8.0;
        Dual x = springSettle(variable(k), m, c, x0, steps, dt);
        real gradAD = x.d;
        real h = 1e-6;
        real fp = springSettle<real>(k + h, m, c, x0, steps, dt);
        real fm = springSettle<real>(k - h, m, c, x0, steps, dt);
        real gradFD = (fp - fm) / (2 * h);
        CHECK_NEAR(gradAD, gradFD, 1e-4);
        CHECK(std::fabs(gradAD) > 1e-6);               // the sim really depends on k
    }

    // ============================ B) ISLAND-PARALLEL SOLVER =====================
    // Many independent islands, each a tiny damped-spring solve. Parallel output
    // must be byte-identical to serial, and (given >1 core) must use >1 thread.
    struct SpringIsland { real k, m, c, x0, dt; int steps; real out; };
    auto solve = [](unsigned i, SpringIsland& is) {
        real x = is.x0, v = 0;
        for (int s = 0; s < is.steps; s++) {
            real a = (-is.k * x - is.c * v) / is.m;
            v += a * is.dt; x += v * is.dt;
        }
        is.out = x + (real)i * (real)1e-9;             // vary per island so slots differ
    };

    std::vector<SpringIsland> base;
    for (int i = 0; i < 2000; i++)
        base.push_back({ 5.0 + 0.01 * i, 1.0, 0.2 + 0.0005 * i, 1.0, 0.01, 300, 0.0 });

    std::vector<SpringIsland> ser = base, par = base;
    solveIslandsSerial(ser, solve);
    unsigned used = solveIslands(par, solve);

    bool identical = true;
    for (size_t i = 0; i < ser.size(); i++)
        if (!(ser[i].out == par[i].out)) { identical = false; break; }   // exact bitwise compare
    CHECK(identical);                                  // parallel == serial, byte for byte
    CHECK(used >= 1);

    unsigned hc = std::thread::hardware_concurrency();
    if (hc > 1) CHECK(used > 1);                       // actually spread across threads
    else std::printf("  (single-core host: thread-count assert skipped, hw=%u)\n", hc);
    std::printf("  parallel solver used %u worker thread(s) over %zu islands\n", used, par.size());

    // ============================ C) CONVEX DECOMPOSITION =======================
    // An L-shaped polygon (one reflex corner) → ≥2 convex pieces tiling the input.
    {
        Polygon2 L = { Vector3(0,0,0), Vector3(2,0,0), Vector3(2,1,0),
                       Vector3(1,1,0), Vector3(1,2,0), Vector3(0,2,0) };
        CHECK(!isConvex(L));                           // it IS concave to begin with
        real areaIn = polygonArea(L);                  // = 3

        std::vector<Polygon2> pieces = convexDecompose(L);
        CHECK(pieces.size() >= 2);

        real areaSum = 0;
        for (const Polygon2& piece : pieces) {
            CHECK(crossSignsAgree(piece));             // every piece is convex
            CHECK(isConvex(piece));                    // (header agrees)
            areaSum += polygonArea(piece);
        }
        CHECK_NEAR(areaSum, areaIn, 1e-9);             // pieces tile the original exactly
        std::printf("  L-shape (area %.3f) → %zu convex pieces (area %.3f)\n",
                    areaIn, pieces.size(), areaSum);
    }

    // C2) a two-reflex "staircase" needs more than one cut; still all convex & tiling.
    {
        Polygon2 S = { Vector3(0,0,0), Vector3(3,0,0), Vector3(3,1,0), Vector3(2,1,0),
                       Vector3(2,2,0), Vector3(1,2,0), Vector3(1,3,0), Vector3(0,3,0) };
        real areaIn = polygonArea(S);
        std::vector<Polygon2> pieces = convexDecompose(S);
        CHECK(pieces.size() >= 3);                      // two concave features → ≥3 pieces
        real areaSum = 0;
        for (const Polygon2& piece : pieces) { CHECK(crossSignsAgree(piece)); areaSum += polygonArea(piece); }
        CHECK_NEAR(areaSum, areaIn, 1e-9);
    }

    // C3) a convex square is returned unchanged (single piece).
    {
        Polygon2 sq = { Vector3(0,0,0), Vector3(1,0,0), Vector3(1,1,0), Vector3(0,1,0) };
        std::vector<Polygon2> pieces = convexDecompose(sq);
        CHECK(pieces.size() == 1);
        CHECK(crossSignsAgree(pieces[0]));
    }

    return test::report("advanced");
}
