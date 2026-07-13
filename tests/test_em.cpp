// Electromagnetics: FDTD (Maxwell curl equations on a Yee grid) and the Method of
// Moments (thin-wire dipole EFIE). Each fact is checked against an analytic or
// classic-textbook result.
#include "phys/fdtd.h"
#include "phys/mom.h"
#include "check.h"
#include <cmath>
#include <vector>
using namespace phys::em;

int main() {
    // ---- FDTD: PEC-cavity eigenfrequency vs analytic f=(c/2)√((1/Lx)²+(1/Ly)²) ----
    {
        int nx = 61, ny = 61; double Lx = nx - 1, Ly = ny - 1;
        FdtdTM s(nx, ny, 0.5); s.pec = true;
        for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++)          // seed the (1,1) mode
            s.Ez[s.id(i, j)] = std::sin(M_PI * i / Lx) * std::sin(M_PI * j / Ly);
        double fana = 0.5 * std::sqrt(1.0 / (Lx * Lx) + 1.0 / (Ly * Ly));
        int pi = 20, pj = 17; double prev = s.Ez[s.id(pi, pj)]; std::vector<double> cross;
        for (int n = 0; n < 6000; n++) { s.step(); double v = s.Ez[s.id(pi, pj)];
            if (prev < 0 && v >= 0) cross.push_back(s.t - v * s.dt / (v - prev)); prev = v; }
        double Tmeas = (cross.back() - cross.front()) / (cross.size() - 1);
        CHECK(cross.size() > 20);
        CHECK_NEAR(1.0 / Tmeas, fana, 0.01 * fana);                        // within 1% (measured ~0.01%)
    }

    // ---- FDTD: a wavefront travels at the speed of light c = 1 ----
    {
        int nx = 400, ny = 40; FdtdTM s(nx, ny, 0.5); s.pec = false; s.border = 0;
        double t1 = -1, t2 = -1; int c1 = 100, c2 = 300;
        for (int n = 0; n < 800; n++) {
            double g = std::exp(-std::pow((n - 40.0) / 12, 2));
            for (int j = 0; j < ny; j++) s.setEz(2, j, g);                 // soft line source
            s.step();
            if (t1 < 0 && std::fabs(s.Ez[s.id(c1, ny / 2)]) > 1e-3) t1 = s.t;
            if (t2 < 0 && std::fabs(s.Ez[s.id(c2, ny / 2)]) > 1e-3) t2 = s.t;
        }
        CHECK(t1 > 0 && t2 > t1);
        CHECK_NEAR((c2 - c1) / (t2 - t1), 1.0, 0.02);                      // < 2% numerical dispersion
    }

    // ---- FDTD: leapfrog is lossless & stable — energy stays bounded in a PEC box ----
    {
        FdtdTM s(60, 60, 0.5); s.pec = true;
        for (int i = 1; i < 59; i++) for (int j = 1; j < 59; j++) {
            double x = (i - 30.0) / 8, y = (j - 30.0) / 8; s.Ez[s.id(i, j)] = std::exp(-(x * x + y * y)); }
        double e0 = s.energy(), emin = e0, emax = e0;
        for (int n = 0; n < 3000; n++) { s.step(); double e = s.energy(); emin = std::min(emin, e); emax = std::max(emax, e); }
        CHECK(emax < 1.10 * e0);                                          // no blow-up
        CHECK(emin > 0.90 * e0);                                          // no numerical dissipation/drift
    }
    // ---- FDTD (TE, the antenna polarization) is likewise stable ----
    {
        FdtdTE s(60, 60, 0.5); s.border = 0;
        for (int i = 1; i < 59; i++) for (int j = 1; j < 59; j++) {
            double x = (i - 30.0) / 8, y = (j - 30.0) / 8; s.Hz[s.id(i, j)] = std::exp(-(x * x + y * y)); }
        double e0 = s.energy(), emin = e0, emax = e0;
        for (int n = 0; n < 3000; n++) { s.step(); double e = s.energy(); emin = std::min(emin, e); emax = std::max(emax, e); }
        CHECK(emax < 1.10 * e0 && emin > 0.90 * e0);
    }

    // ---- MoM: near-resonant half-wave dipole matches the classic ~73 Ω ----
    {
        Dipole d; d.L = 0.46; d.a = 0.005; d.N = 41; d.solve();
        CHECK(d.Zin.real() > 65 && d.Zin.real() < 85);                    // R ≈ 73 Ω
        CHECK(std::fabs(d.Zin.imag()) < 12);                              // near resonance → small reactance
        // current distribution is symmetric and peaks at the feed
        double asym = 0; for (int k = 0; k < d.N; k++) asym = std::max(asym, std::fabs(std::abs(d.I[k]) - std::abs(d.I[d.N - 1 - k])));
        CHECK(asym < 1e-9);
        int am = 0; for (int k = 0; k < d.N; k++) if (std::abs(d.I[k]) > std::abs(d.I[am])) am = k;
        CHECK(std::abs(am - d.N / 2) <= 2);                               // maximum at (or beside) the feed
        CHECK(std::abs(d.I[0]) < 0.3 * d.maxCurrent());                   // current tapers toward the ends
    }
    // ---- MoM: reactance is capacitive when short, inductive at 0.5λ (classic behavior) ----
    {
        Dipole sht; sht.L = 0.40; sht.a = 0.005; sht.N = 41; sht.solve();
        Dipole hlf; hlf.L = 0.50; hlf.a = 0.005; hlf.N = 41; hlf.solve();
        CHECK(sht.Zin.imag() < 0);                                        // short dipole: capacitive
        CHECK(hlf.Zin.imag() > 0);                                        // full 0.5λ: inductive
        CHECK(hlf.Zin.real() > 0);                                        // passive
    }
    // ---- MoM: radiation pattern has a null along the wire axis and a broadside max ----
    {
        Dipole d; d.L = 0.5; d.a = 0.005; d.N = 41; d.solve();
        CHECK(d.pattern(0.0) < 1e-6 * d.pattern(M_PI / 2));               // deep null at θ=0
        CHECK(d.pattern(M_PI / 2) > d.pattern(M_PI / 4));                 // broadside is the maximum
        CHECK(d.pattern(M_PI / 4) > 0.4 * d.pattern(M_PI / 2));           // smooth, broad main lobe
    }

    return test::report("em");
}
