// Non-Newtonian fluids via Smoothed-Particle Hydrodynamics (SPH).
//
// A weakly-compressible SPH fluid (Müller et al. 2003 kernels: poly6 density,
// spiky pressure gradient, viscosity Laplacian) whose viscosity is *shear-rate
// dependent* — a generalised-Newtonian (Herschel–Bulkley / power-law) rheology
//     μ(γ̇) = μ_min + K·γ̇^(n-1) + τ_y·(1 − e^(−m γ̇))/γ̇          (Papanastasiou-regularised yield)
// with the local shear rate γ̇ = √(2 D:D) taken from the SPH velocity-gradient
// strain-rate tensor D. This one knob reproduces the whole zoo:
//   • n > 1  → shear-thickening (dilatant): oobleck — solid under fast impact, flows when slow
//   • n < 1  → shear-thinning  (pseudoplastic): ketchup, paint, blood
//   • τ_y>0  → yield stress (Bingham/Herschel–Bulkley): toothpaste, mayonnaise
//   • n = 1, τ_y = 0 → ordinary Newtonian fluid (water)
// A rigid sphere can be dropped in and is two-way coupled, so a thickening fluid
// visibly resists a fast impactor while a Newtonian one lets it plunge.
#pragma once
#include "core.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace phys {

struct Rheology {
    real muMin = 0.1;      // base (infinite-shear) viscosity
    real K = 3.5;          // consistency index
    real n = 1.0;          // flow index: >1 thickening, <1 thinning, =1 Newtonian
    real tauY = 0.0;       // yield stress (Bingham/Herschel–Bulkley)
    real mReg = 200.0;     // yield regularisation sharpness
    real muMax = 60.0;     // clamp so the explicit solver stays stable

    real viscosity(real gammaDot) const {
        real g = gammaDot;
        real mu = muMin + K * real_pow(g > 1e-6 ? g : 1e-6, n - 1);
        if (tauY > 0) mu += tauY * (1 - std::exp(-mReg * g)) / (g + 1e-6);
        return mu < muMin ? muMin : (mu > muMax ? muMax : mu);
    }

    static Rheology water()    { return {0.5, 3.5, 1.0, 0.0, 200.0, 60.0}; }
    static Rheology oobleck()  { return {0.4, 9.0, 2.2, 0.0, 200.0, 900.0}; }   // shear-thickening
    static Rheology ketchup()  { return {0.3, 6.0, 0.4, 12.0, 250.0, 200.0}; }  // shear-thinning + yield
    static Rheology lava()     { return {6.0, 30.0, 0.7, 9.0, 60.0, 400.0}; }   // hot Bingham: thick, shear-thinning, yields
};

// A fixed vertical cylinder obstacle (axis +y) the fluid must flow around.
struct SPHCylinder { real cx = 0, cz = 0, radius = 0.15, y0 = 0, y1 = 1e9, friction = 0.8; };

// A rigid sphere that floats/sinks in the fluid (two-way coupled).
struct FluidBall {
    Vector3 pos, vel; real radius = 0.1, mass = 1.0; bool active = false;
    Vector3 force;                                   // reaction accumulated from the fluid each step
};

class SPHFluid {
public:
    real h = 0.10;                                   // smoothing radius
    real mass = 0.125;                               // per-particle mass
    real rho0 = 1000;                                // rest density (measured in build())
    real stiffness = 200;                            // EOS: p = stiffness·(ρ − ρ0)
    Vector3 gravity = Vector3(0, -9.81, 0);
    real restitution = 0.3, wallFriction = 0.5;      // container walls
    real xsphEps = 0.18;                             // XSPH velocity smoothing (kills SPH jitter)
    Rheology rheo;

    Vector3 boundsMin, boundsMax;                    // container AABB
    std::vector<Vector3> pos, vel;
    std::vector<real> rho, pres, gammaDot;
    std::vector<FluidBall> balls;
    std::vector<SPHCylinder> cylinders;

    // Add one particle at runtime (for inflow/emitters); keeps all arrays in sync.
    void addParticle(const Vector3& p, const Vector3& v) {
        pos.push_back(p); vel.push_back(v);
        rho.push_back(rho0); pres.push_back(0); gammaDot.push_back(0);
    }
    // Swap-remove particle i (for outflow/drains). Any parallel per-particle array
    // the caller keeps must mirror this: arr[i] = arr.back(); arr.pop_back();
    void removeParticle(int i) {
        int last = (int)pos.size() - 1;
        if (i != last) { pos[i] = pos[last]; vel[i] = vel[last]; rho[i] = rho[last]; pres[i] = pres[last]; gammaDot[i] = gammaDot[last]; }
        pos.pop_back(); vel.pop_back(); rho.pop_back(); pres.pop_back(); gammaDot.pop_back();
    }

    // --- kernels (3D) ---
    real poly6(real r2) const {                      // pass squared distance
        if (r2 >= h2) return 0; real d = h2 - r2; return cPoly6 * d * d * d;
    }
    Vector3 spikyGrad(const Vector3& rij, real r) const {   // ∇W_spiky(r_i − r_j)
        if (r >= h || r < 1e-9) return Vector3();
        real c = cSpiky * (h - r) * (h - r) / r; return rij * c;      // points along rij
    }
    real viscLap(real r) const { return r < h ? cVisc * (h - r) : 0; }

    void setKernelConstants() {
        h2 = h * h;
        cPoly6 = 315.0 / (64.0 * real_pi * real_pow(h, 9));
        cSpiky = -45.0 / (real_pi * real_pow(h, 6));
        cVisc = 45.0 / (real_pi * real_pow(h, 6));
    }

    // Fill an axis-aligned block of fluid inside a container; measures rho0 so the
    // resting fluid starts at ~zero pressure regardless of spacing/mass.
    void build(const Vector3& cmin, const Vector3& cmax,
               const Vector3& fmin, const Vector3& fmax, real spacing) {
        boundsMin = cmin; boundsMax = cmax; setKernelConstants();
        pos.clear(); vel.clear();
        for (real x = fmin.x; x <= fmax.x; x += spacing)
            for (real y = fmin.y; y <= fmax.y; y += spacing)
                for (real z = fmin.z; z <= fmax.z; z += spacing) {
                    pos.push_back(Vector3(x, y, z)); vel.push_back(Vector3());
                }
        rho.assign(pos.size(), rho0); pres.assign(pos.size(), 0); gammaDot.assign(pos.size(), 0);
        buildGrid(); computeDensity();
        real sum = 0; int cnt = 0;                   // rest density = mean interior density
        for (size_t i = 0; i < pos.size(); i++) if (rho[i] > 0) { sum += rho[i]; cnt++; }
        if (cnt) rho0 = sum / cnt;
    }

    int nParticles() const { return (int)pos.size(); }

    void step(real dt) {
        curDt = dt;
        buildGrid();
        computeDensity();
        for (size_t i = 0; i < pos.size(); i++) {
            real p = stiffness * (rho[i] - rho0); pres[i] = p > 0 ? p : 0;   // clamp ≥0 (no tensile clumping)
        }
        computeShearRate();
        std::vector<Vector3> force(pos.size());
        computeForces(force);
        // integrate (semi-implicit Euler) + boundaries + ball coupling
        for (auto& b : balls) if (b.active) b.force = Vector3();
        for (size_t i = 0; i < pos.size(); i++) {
            Vector3 acc = force[i] * (1.0 / rho[i]) + gravity;
            vel[i] += acc * dt;
        }
        // XSPH: advect with a neighbour-smoothed velocity (momentum-preserving)
        std::vector<Vector3> vadv(pos.size());
        for (int i = 0; i < (int)pos.size(); i++) {
            Vector3 corr, pi = pos[i], vi = vel[i];
            forNeighbours(pi, [&](int j) {
                if (j == i) return; real r2 = (pi - pos[j]).squareMagnitude();
                if (r2 < h2) corr += (vel[j] - vi) * (mass / rho[j] * poly6(r2));
            });
            vadv[i] = vi + corr * xsphEps;
        }
        for (size_t i = 0; i < pos.size(); i++) {
            pos[i] += vadv[i] * dt;
            resolveBalls((int)i);
            resolveCylinders((int)i);
            resolveWalls((int)i);
        }
        for (auto& b : balls) if (b.active) {
            b.vel += (b.force * (1.0 / b.mass) + gravity) * dt;
            b.pos += b.vel * dt;
            // ball vs container floor/walls
            for (int a = 0; a < 3; a++) {
                if (b.pos[a] - b.radius < boundsMin[a]) { setAxis(b.pos, a, boundsMin[a] + b.radius); flip(b.vel, a, 0.2); }
                if (b.pos[a] + b.radius > boundsMax[a]) { setAxis(b.pos, a, boundsMax[a] - b.radius); flip(b.vel, a, 0.2); }
            }
        }
    }

    real kineticEnergy() const {
        real e = 0; for (size_t i = 0; i < pos.size(); i++) e += 0.5 * mass * vel[i].squareMagnitude(); return e;
    }

private:
    real curDt = 1;
    real h2 = 0.01, cPoly6 = 1, cSpiky = 1, cVisc = 1;
    // uniform grid neighbour search over the container
    int gnx = 1, gny = 1, gnz = 1; Vector3 gorigin; real gcell = 0.1;
    std::vector<std::vector<int>> cells;

    void buildGrid() {
        gcell = h; gorigin = boundsMin;
        gnx = std::max(1, (int)std::ceil((boundsMax.x - boundsMin.x) / gcell) + 1);
        gny = std::max(1, (int)std::ceil((boundsMax.y - boundsMin.y) / gcell) + 1);
        gnz = std::max(1, (int)std::ceil((boundsMax.z - boundsMin.z) / gcell) + 1);
        cells.assign((size_t)gnx * gny * gnz, {});
        for (int i = 0; i < (int)pos.size(); i++) { int c = cellIndex(pos[i]); if (c >= 0) cells[c].push_back(i); }
    }
    int clampi(int v, int hi) const { return v < 0 ? 0 : (v > hi ? hi : v); }
    int cellIndex(const Vector3& p) const {
        int cx = clampi((int)((p.x - gorigin.x) / gcell), gnx - 1);
        int cy = clampi((int)((p.y - gorigin.y) / gcell), gny - 1);
        int cz = clampi((int)((p.z - gorigin.z) / gcell), gnz - 1);
        return (cx * gny + cy) * gnz + cz;
    }
    // call fn(j) for every particle j within h of position p
    template <class F> void forNeighbours(const Vector3& p, F&& fn) const {
        int cx = clampi((int)((p.x - gorigin.x) / gcell), gnx - 1);
        int cy = clampi((int)((p.y - gorigin.y) / gcell), gny - 1);
        int cz = clampi((int)((p.z - gorigin.z) / gcell), gnz - 1);
        for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
            int x = cx + dx, y = cy + dy, z = cz + dz;
            if (x < 0 || y < 0 || z < 0 || x >= gnx || y >= gny || z >= gnz) continue;
            for (int j : cells[(x * gny + y) * gnz + z]) fn(j);
        }
    }

    void computeDensity() {
        for (int i = 0; i < (int)pos.size(); i++) {
            real d = 0; Vector3 pi = pos[i];
            forNeighbours(pi, [&](int j) { real r2 = (pi - pos[j]).squareMagnitude(); if (r2 < h2) d += mass * poly6(r2); });
            rho[i] = d > 1e-6 ? d : 1e-6;
        }
    }

    // local shear rate from the SPH velocity-gradient strain-rate tensor
    void computeShearRate() {
        for (int i = 0; i < (int)pos.size(); i++) {
            real g[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            Vector3 pi = pos[i], vi = vel[i];
            forNeighbours(pi, [&](int j) {
                if (j == i) return; Vector3 rij = pi - pos[j]; real r = rij.magnitude();
                if (r >= h || r < 1e-9) return;
                Vector3 gw = spikyGrad(rij, r);              // ∇_i W_ij
                Vector3 dv = vel[j] - vi;                    // (v_j − v_i)
                real w = mass / rho[j];
                for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) g[3 * a + b] += w * dv[a] * gw[b];
            });
            real gd = 0;                                     // γ̇ = sqrt(2 D:D), D = ½(∇v+∇vᵀ)
            for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
                real D = 0.5 * (g[3 * a + b] + g[3 * b + a]); gd += 2 * D * D;
            }
            gammaDot[i] = real_sqrt(gd);
        }
    }

    void computeForces(std::vector<Vector3>& force) {
        std::vector<real> mu(pos.size());
        for (size_t i = 0; i < pos.size(); i++) mu[i] = rheo.viscosity(gammaDot[i]);
        for (int i = 0; i < (int)pos.size(); i++) {
            Vector3 fp, fv; Vector3 pi = pos[i], vi = vel[i]; real presi = pres[i];
            forNeighbours(pi, [&](int j) {
                if (j == i) return; Vector3 rij = pi - pos[j]; real r = rij.magnitude();
                if (r >= h || r < 1e-9) return;
                Vector3 gw = spikyGrad(rij, r);
                fp += gw * (-mass * (presi + pres[j]) / (2 * rho[j]));         // pressure
                real muij = 0.5 * (mu[i] + mu[j]);
                fv += (vel[j] - vi) * (muij * mass / rho[j] * viscLap(r));     // shear-dependent viscosity
            });
            force[i] = fp + fv;
        }
    }

    // push particle i out of any ball and exchange normal momentum (two-way)
    void resolveBalls(int i) {
        for (auto& b : balls) {
            if (!b.active) continue;
            Vector3 d = pos[i] - b.pos; real dist = d.magnitude();
            if (dist < b.radius && dist > 1e-9) {
                Vector3 nrm = d * (1.0 / dist);
                pos[i] = b.pos + nrm * b.radius;
                real vrel = (vel[i] - b.vel) * nrm;
                if (vrel < 0) {
                    Vector3 dp = nrm * (vrel * mass);        // momentum removed from the particle (points inward)
                    vel[i] -= dp * (1.0 / mass);             // particle bounces off the ball surface
                    b.force += dp * (1.0 / curDt);           // equal-and-opposite reaction lifts the ball
                }
            }
        }
    }
    void resolveCylinders(int i) {
        for (auto& c : cylinders) {
            if (pos[i].y < c.y0 || pos[i].y > c.y1) continue;
            real dx = pos[i].x - c.cx, dz = pos[i].z - c.cz, d = real_sqrt(dx * dx + dz * dz);
            if (d < c.radius && d > 1e-9) {
                real nx = dx / d, nz = dz / d;
                pos[i].x = c.cx + nx * c.radius; pos[i].z = c.cz + nz * c.radius;
                real vn = vel[i].x * nx + vel[i].z * nz;
                if (vn < 0) { vel[i].x -= vn * nx; vel[i].z -= vn * nz;   // kill inward normal velocity
                    vel[i].x *= c.friction; vel[i].z *= c.friction; }
            }
        }
    }
    void resolveWalls(int i) {
        for (int a = 0; a < 3; a++) {
            if (pos[i][a] < boundsMin[a]) { setAxis(pos[i], a, boundsMin[a]); flip(vel[i], a, restitution); dampTangent(vel[i], a); }
            if (pos[i][a] > boundsMax[a]) { setAxis(pos[i], a, boundsMax[a]); flip(vel[i], a, restitution); dampTangent(vel[i], a); }
        }
    }
    static void setAxis(Vector3& v, int a, real val) { if (a == 0) v.x = val; else if (a == 1) v.y = val; else v.z = val; }
    static void flip(Vector3& v, int a, real rest) { real c = v[a]; setAxis(v, a, -c * rest); }
    void dampTangent(Vector3& v, int a) { for (int b = 0; b < 3; b++) if (b != a) setAxis(v, b, v[b] * wallFriction); }
};

} // namespace phys
