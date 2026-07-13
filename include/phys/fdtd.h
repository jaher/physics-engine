// FDTD — Finite-Difference Time-Domain solution of Maxwell's curl equations on a
// Yee grid (Yee 1966; Taflove). The electric and magnetic fields are staggered in
// both space and time and leapfrogged forward, which reproduces wave propagation,
// interference, and radiation directly from Ampère's and Faraday's laws with no
// analytic ansatz. Two 2-D polarizations are provided:
//
//   FdtdTM  (transverse-magnetic to z):  Ez, Hx, Hy   — clean PEC-cavity modes,
//           used to validate wave speed / resonance / energy conservation.
//   FdtdTE  (transverse-electric to z):  Ex, Ey, Hz   — the in-plane E field, used
//           to drive and render a dipole antenna (current Jy → radiated E,H).
//
// Units are normalized (c = ε₀ = μ₀ = 1, Δx = Δy = 1); the Courant number Sc < 1/√2
// sets the timestep dt = Sc. Boundaries are either PEC (fields vanish → a lossless
// resonant box) or a graded absorbing border that soaks outgoing waves for
// open-region radiation problems.
#pragma once
#include <vector>
#include <cmath>

namespace phys { namespace em {

// ---- transverse-magnetic (Ez, Hx, Hy) --------------------------------------
struct FdtdTM {
    int nx, ny; double dt;
    bool pec = true;              // PEC walls (lossless box) vs absorbing border
    int  border = 0;              // absorbing-layer width when !pec
    std::vector<double> Ez, Hx, Hy;
    double t = 0;

    FdtdTM(int nx_, int ny_, double Sc = 0.5) : nx(nx_), ny(ny_), dt(Sc) {
        Ez.assign(nx * ny, 0.0); Hx.assign(nx * ny, 0.0); Hy.assign(nx * ny, 0.0);
    }
    int id(int i, int j) const { return i * ny + j; }

    void step() {
        // Faraday: H at n+1/2
        for (int i = 0; i < nx; i++) for (int j = 0; j < ny - 1; j++)
            Hx[id(i, j)] -= dt * (Ez[id(i, j + 1)] - Ez[id(i, j)]);           // ∂Hx/∂t = -∂Ez/∂y
        for (int i = 0; i < nx - 1; i++) for (int j = 0; j < ny; j++)
            Hy[id(i, j)] += dt * (Ez[id(i + 1, j)] - Ez[id(i, j)]);           // ∂Hy/∂t = +∂Ez/∂x
        // Ampère: Ez at n+1 (interior; boundary handled below)
        for (int i = 1; i < nx - 1; i++) for (int j = 1; j < ny - 1; j++)
            Ez[id(i, j)] += dt * ((Hy[id(i, j)] - Hy[id(i - 1, j)]) - (Hx[id(i, j)] - Hx[id(i, j - 1)]));
        if (!pec) absorb();
        t += dt;
    }
    // graded absorbing border: bleed the fields toward zero within `border` cells
    void absorb() {
        if (border <= 0) return;
        for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) {
            int d = std::min(std::min(i, nx - 1 - i), std::min(j, ny - 1 - j));
            if (d < border) { double s = 1.0 - 0.5 * (double)(border - d) / border * (double)(border - d) / border;
                Ez[id(i, j)] *= s; Hx[id(i, j)] *= s; Hy[id(i, j)] *= s; }
        }
    }
    void addEz(int i, int j, double v) { Ez[id(i, j)] += v; }   // soft source
    void setEz(int i, int j, double v) { Ez[id(i, j)] = v; }    // hard source
    double energy() const {                                     // ½∫(E²+H²)
        double e = 0; for (int k = 0; k < nx * ny; k++) e += Ez[k] * Ez[k] + Hx[k] * Hx[k] + Hy[k] * Hy[k];
        return 0.5 * e;
    }
};

// ---- transverse-electric (Ex, Ey, Hz) --------------------------------------
struct FdtdTE {
    int nx, ny; double dt;
    int border = 12;                        // absorbing border (open region by default)
    std::vector<double> Ex, Ey, Hz;
    double t = 0;

    FdtdTE(int nx_, int ny_, double Sc = 0.5) : nx(nx_), ny(ny_), dt(Sc) {
        Ex.assign(nx * ny, 0.0); Ey.assign(nx * ny, 0.0); Hz.assign(nx * ny, 0.0);
    }
    int id(int i, int j) const { return i * ny + j; }

    void updateH() {
        for (int i = 0; i < nx - 1; i++) for (int j = 0; j < ny - 1; j++)
            Hz[id(i, j)] += dt * ((Ex[id(i, j + 1)] - Ex[id(i, j)]) - (Ey[id(i + 1, j)] - Ey[id(i, j)]));
    }
    void updateE() {
        for (int i = 0; i < nx; i++) for (int j = 1; j < ny; j++)
            Ex[id(i, j)] += dt * (Hz[id(i, j)] - Hz[id(i, j - 1)]);           // ∂Ex/∂t = +∂Hz/∂y
        for (int i = 1; i < nx; i++) for (int j = 0; j < ny; j++)
            Ey[id(i, j)] -= dt * (Hz[id(i, j)] - Hz[id(i - 1, j)]);           // ∂Ey/∂t = -∂Hz/∂x
    }
    void absorb() {
        if (border <= 0) return;
        for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) {
            int d = std::min(std::min(i, nx - 1 - i), std::min(j, ny - 1 - j));
            if (d < border) { double r = (double)(border - d) / border; double s = 1.0 - 0.35 * r * r;
                Ex[id(i, j)] *= s; Ey[id(i, j)] *= s; Hz[id(i, j)] *= s; }
        }
    }
    // one full leapfrog step; injectJy(cellId)→current density added to Ey (the feed)
    template <class Src>
    void step(Src&& injectJy) {
        updateH();
        updateE();
        for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) {
            double J = injectJy(i, j); if (J != 0.0) Ey[id(i, j)] -= dt * J;
        }
        absorb();
        t += dt;
    }
    void step() { step([](int, int) { return 0.0; }); }
    double energy() const {
        double e = 0; for (int k = 0; k < nx * ny; k++) e += Ex[k] * Ex[k] + Ey[k] * Ey[k] + Hz[k] * Hz[k];
        return 0.5 * e;
    }
};

}} // namespace phys::em
