// Spring harmonics: a coupled spring-mass chain must ring at the analytic
// normal-mode frequencies ω_n = 2√(k/m)·sin(nπ/2(P+1)), and conserve energy.
#include "phys/springs.h"
#include "check.h"
#include <vector>
#include <algorithm>
using namespace phys;

// Measure the oscillation period of an interior node by timing upward zero
// crossings of its longitudinal displacement.
static real measurePeriod(SpringChain& c, int watch, real dt, int steps) {
    std::vector<real> crossings;
    real prev = c.displacement(watch).x;
    for (int i = 0; i < steps; i++) {
        c.step(dt);
        real d = c.displacement(watch).x;
        if (prev <= 0 && d > 0) crossings.push_back(i * dt);       // rising through zero
        prev = d;
    }
    if (crossings.size() < 2) return 0;
    return (crossings.back() - crossings.front()) / (crossings.size() - 1);
}

int main() {
    // A) Longitudinal mode frequencies match the analytic dispersion relation.
    //    (Longitudinal motion is exactly linear, so this is a sharp test.)
    {
        for (int n = 1; n <= 3; n++) {
            SpringChain c;
            c.build(/*P*/5, /*spacing*/1.0, /*k*/100.0, /*m*/1.0, Vector3(0, 0, 0));
            c.setMode(n, 0.02, c.axis);                            // longitudinal pluck
            real Tsim = measurePeriod(c, 1, 1e-4, 200000);
            real Tref = c.modePeriod(n);
            CHECK(Tsim > 0);
            CHECK_NEAR(Tsim, Tref, Tref * 0.01);                   // within 1%
        }
    }

    // B) Single interior mass between two walls: ω = √(2k/m) (two springs in parallel).
    {
        SpringChain c;
        c.build(1, 1.0, 100.0, 1.0, Vector3(0, 0, 0));
        CHECK_NEAR(c.modeOmega(1), real_sqrt(2 * 100.0 / 1.0), 1e-9);
    }

    // C) Modes are ordered and the fundamental is the slowest (harmonic series).
    {
        SpringChain c;
        c.build(8, 1.0, 50.0, 1.0, Vector3(0, 0, 0));
        for (int n = 1; n < 8; n++) CHECK(c.modeOmega(n + 1) > c.modeOmega(n));
    }

    // D) Energy is conserved (undamped): a ringing chain holds its total energy.
    {
        SpringChain c;
        c.build(8, 1.0, 80.0, 1.0, Vector3(0, 0, 0));
        c.setMode(2, 0.05, c.axis);
        real E0 = c.energy(), lo = E0, hi = E0;
        for (int i = 0; i < 40000; i++) { c.step(5e-4); real E = c.energy(); lo = std::min(lo, E); hi = std::max(hi, E); }
        CHECK(E0 > 0);
        CHECK((hi - lo) / E0 < 0.02);                              // < 2% drift over ~many periods
    }

    // E) Transverse harmonics: a pre-tensioned chain gives a near-linear harmonic
    //    series (ω_n ≈ n·ω_1) — the plucked-string spectrum.
    {
        SpringChain c;
        c.build(40, 0.1, 200.0, 0.02, Vector3(0, 0, 0), Vector3(1, 0, 0), /*restLength*/0.05);
        real w1 = c.modeOmega(1, true);
        for (int n : {2, 3, 4}) CHECK_NEAR(c.modeOmega(n, true) / w1, (real)n, 0.05 * n);
    }

    return test::report("harmonics");
}
