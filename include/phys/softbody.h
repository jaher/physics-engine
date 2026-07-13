// Volumetric deformables and position-based fluids — the two families the engine
// was missing (FEATURES.md listed "soft-body volumes/FEM" and a PBF fluid as
// not-implemented). All header-only, built on the engine's Vector3 / Matrix3.
//
//   1. SoftBody   — a tetrahedral CO-ROTATIONAL linear-FEM deformable
//                   (Müller & Gross, "Interactive Virtual Materials", 2004).
//                   Per element the deformation gradient F is factored F = R·S by
//                   polar decomposition; the linear-elastic internal force is
//                   evaluated in the un-rotated frame and rotated back —
//                   f = −R·K·(Rᵀx − x0) — which removes the ghost forces a plain
//                   linear model produces under large rotations. Builds a solid
//                   box, drops it on the ground, and it squashes and rebounds
//                   while conserving volume.
//   2. resolveSelfCollisions — a uniform-spatial-hash particle de-overlap helper
//                   so phys::Cloth (or any point set) can self-collide without
//                   touching cloth.h.
//   3. PBFFluid   — Position-Based Fluids (Macklin & Müller, 2013): poly6/spiky
//                   kernels, a per-particle density-constraint multiplier λ and
//                   iterative position projection for incompressibility. A
//                   projective, unconditionally-stable alternative to the SPH in
//                   sph.h.
#pragma once
#include "core.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <unordered_map>

namespace phys {

// ------------------------------------------------------------------ shared math
// 3x3 determinant (Matrix3 exposes an inverse but not a determinant getter).
inline real det3(const Matrix3& m) {
    const real* d = m.data;
    return d[0] * (d[4] * d[8] - d[5] * d[7])
         - d[1] * (d[3] * d[8] - d[5] * d[6])
         + d[2] * (d[3] * d[7] - d[4] * d[6]);
}

// Orthogonal (rotation) factor R of the polar decomposition F = R·S, via the
// quadratically-convergent Newton/Higham iteration  X ← ½(X + X⁻ᵀ),  X₀ = F.
//   • F a pure rotation  → converges in one step, returns F unchanged.
//   • F symmetric (pure stretch) → R → I (all the "rotation" is really stretch).
inline Matrix3 polarRotation(const Matrix3& F, int maxIter = 64) {
    Matrix3 X = F;
    for (int it = 0; it < maxIter; ++it) {
        if (real_abs(det3(X)) < (real)1e-12) break;         // singular: bail out
        Matrix3 Xit = X.inverse().transpose();              // X⁻ᵀ
        Matrix3 Xn = X; Xn += Xit; Xn *= (real)0.5;
        real diff = 0;
        for (int n = 0; n < 9; ++n) { real e = Xn.data[n] - X.data[n]; diff += e * e; }
        X = Xn;
        if (diff < (real)1e-24) break;
    }
    return X;
}
// Full polar decomposition F = R·S (R orthogonal, S symmetric positive-definite).
inline void polarDecompose(const Matrix3& F, Matrix3& R, Matrix3& S) {
    R = polarRotation(F);
    S = R.transpose() * F;
}

// ============================================================= 1. SoftBody (FEM)
// A tetrahedral co-rotational linear-FEM deformable. Mass is lumped to the nodes;
// each tetrahedron caches its rest inverse-shape matrix Dm⁻¹ and its constant
// 12×12 element stiffness Kₑ = V·Bᵀ·C·B, assembled from the isotropic Lamé
// constants of the material.
class SoftBody {
public:
    struct Node { Vector3 x, v; real mass = 0, invMass = 0; };
    struct Tet {
        int n[4];                 // node indices
        Matrix3 DmInv;            // inverse rest edge matrix
        real Ke[12][12];          // element stiffness (constant, un-rotated frame)
        real restVol = 0;
    };

    std::vector<Node> nodes;
    std::vector<Vector3> rest;    // rest (material) position of each node
    std::vector<Tet> tets;

    Vector3 gravity = Vector3(0, -9.81, 0);
    real groundY     = 0;
    real restitution = (real)0.3;   // ground bounce
    real friction    = (real)0.3;   // tangential grip at the ground
    real damping     = (real)0.999; // global velocity damping (explicit-solver sink)

    // element stiffness block  K_ij = V(λ bᵢbⱼᵀ + μ(bᵢ·bⱼ)I + μ bⱼbᵢᵀ)
    void addTet(int a, int b, int c, int d, real lambda, real mu, real density) {
        Vector3 p0 = rest[a], p1 = rest[b], p2 = rest[c], p3 = rest[d];
        Matrix3 Dm; Dm.setComponents(p1 - p0, p2 - p0, p3 - p0);
        real vol = det3(Dm) / 6;
        if (vol < 0) {                                  // keep positive orientation
            std::swap(c, d); std::swap(p2, p3);
            Dm.setComponents(p1 - p0, p2 - p0, p3 - p0);
            vol = det3(Dm) / 6;
        }
        if (vol < (real)1e-12) return;                  // degenerate tet
        Tet t; t.n[0] = a; t.n[1] = b; t.n[2] = c; t.n[3] = d;
        t.DmInv = Dm.inverse(); t.restVol = vol;
        // shape-function gradients: rows of Dm⁻¹ are ∇λ₁,∇λ₂,∇λ₃; ∇λ₀ = −Σ.
        Vector3 r0 = t.DmInv.getRowVector(0), r1 = t.DmInv.getRowVector(1), r2 = t.DmInv.getRowVector(2);
        Vector3 B[4] = { -(r0 + r1 + r2), r0, r1, r2 };
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                real dot = B[i] * B[j];
                for (int x = 0; x < 3; ++x)
                    for (int y = 0; y < 3; ++y) {
                        real val = lambda * B[i][x] * B[j][y] + mu * B[j][x] * B[i][y];
                        if (x == y) val += mu * dot;
                        t.Ke[3 * i + x][3 * j + y] = vol * val;
                    }
            }
        tets.push_back(t);
        real m = density * vol / 4;                     // lumped mass
        nodes[a].mass += m; nodes[b].mass += m; nodes[c].mass += m; nodes[d].mass += m;
    }

    void finalizeMasses() {
        for (auto& n : nodes) n.invMass = n.mass > 0 ? ((real)1) / n.mass : 0;
    }
    void pin(int i) { nodes[i].invMass = 0; }

    // Solid box of nx×ny×nz cells (each a cube of edge `size`), split into 6 tets
    // per cell (Kuhn triangulation about the 0–7 diagonal so faces match between
    // cells). Young's modulus / Poisson ratio give the Lamé constants.
    static SoftBody makeBox(int nx, int ny, int nz, real size,
                            real youngs, real poisson, real density,
                            const Vector3& origin = Vector3(0, 0, 0)) {
        SoftBody sb;
        int NX = nx + 1, NY = ny + 1, NZ = nz + 1;
        auto id = [=](int i, int j, int k) { return (i * NY + j) * NZ + k; };
        sb.nodes.resize(NX * NY * NZ);
        sb.rest.resize(NX * NY * NZ);
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    Vector3 p = origin + Vector3(i * size, j * size, k * size);
                    sb.nodes[id(i, j, k)].x = p;
                    sb.rest[id(i, j, k)]    = p;
                }
        real lambda = youngs * poisson / ((1 + poisson) * (1 - 2 * poisson));
        real mu     = youngs / (2 * (1 + poisson));
        static const int off[8][3] = {
            {0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
        static const int cellTet[6][4] = {
            {0,7,1,3},{0,7,3,2},{0,7,2,6},{0,7,6,4},{0,7,4,5},{0,7,5,1}};
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                for (int k = 0; k < nz; ++k) {
                    int corner[8];
                    for (int c = 0; c < 8; ++c)
                        corner[c] = id(i + off[c][0], j + off[c][1], k + off[c][2]);
                    for (int t = 0; t < 6; ++t)
                        sb.addTet(corner[cellTet[t][0]], corner[cellTet[t][1]],
                                  corner[cellTet[t][2]], corner[cellTet[t][3]],
                                  lambda, mu, density);
                }
        sb.finalizeMasses();
        return sb;
    }

    // One semi-implicit Euler step of the co-rotational FEM.
    void step(real dt) {
        std::vector<Vector3> force(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) force[i] = gravity * nodes[i].mass;
        for (const Tet& t : tets) {
            Vector3 x0 = nodes[t.n[0]].x, x1 = nodes[t.n[1]].x,
                    x2 = nodes[t.n[2]].x, x3 = nodes[t.n[3]].x;
            Matrix3 Ds; Ds.setComponents(x1 - x0, x2 - x0, x3 - x0);
            Matrix3 F = Ds * t.DmInv;                    // deformation gradient
            Matrix3 R = polarRotation(F);                // extract rotation
            Vector3 d[4];                                // un-rotated displacement Rᵀx − x0
            for (int i = 0; i < 4; ++i)
                d[i] = R.transformTranspose(nodes[t.n[i]].x) - rest[t.n[i]];
            for (int i = 0; i < 4; ++i) {                // y_i = Σ_j Kₑ_ij d_j
                Vector3 yi;
                for (int j = 0; j < 4; ++j)
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b)
                            yi[a] += t.Ke[3 * i + a][3 * j + b] * d[j][b];
                force[t.n[i]] += R.transform(yi) * (real)(-1);   // f_i = −R·y_i (restoring)
            }
        }
        for (size_t i = 0; i < nodes.size(); ++i) {
            Node& nd = nodes[i];
            if (nd.invMass == 0) { nd.v = Vector3(); continue; }
            nd.v += force[i] * (nd.invMass * dt);
            nd.v *= damping;
            nd.x += nd.v * dt;
            if (nd.x.y < groundY) {                      // ground-plane collision
                nd.x.y = groundY;
                if (nd.v.y < 0) nd.v.y = -nd.v.y * restitution;
                nd.v.x *= (1 - friction);
                nd.v.z *= (1 - friction);
            }
        }
    }

    // Current deformed volume = Σ |det(edge matrix)|/6 over all tets.
    real totalVolume() const {
        real V = 0;
        for (const Tet& t : tets) {
            Vector3 x0 = nodes[t.n[0]].x, x1 = nodes[t.n[1]].x,
                    x2 = nodes[t.n[2]].x, x3 = nodes[t.n[3]].x;
            Matrix3 Ds; Ds.setComponents(x1 - x0, x2 - x0, x3 - x0);
            V += real_abs(det3(Ds)) / 6;
        }
        return V;
    }
    real restVolume() const { real V = 0; for (const Tet& t : tets) V += t.restVol; return V; }
    Vector3 centroid() const {
        Vector3 c; for (const Node& n : nodes) c += n.x;
        if (!nodes.empty()) c *= ((real)1) / (real)nodes.size();
        return c;
    }
};

// =================================================== 2. cloth self-collision help
// Push apart every particle pair closer than 2·radius, using a uniform spatial
// hash (cell size = 2·radius, so only the 27 neighbouring cells need checking).
// Symmetric half-corrections; a couple of passes firm it up. Lets phys::Cloth (or
// any point cloud) self-collide without editing cloth.h.
inline void resolveSelfCollisions(std::vector<Vector3>& pos, real radius, int iterations = 2) {
    if (radius <= 0 || pos.size() < 2) return;
    real diam = 2 * radius, cell = diam, inv = ((real)1) / cell;
    auto key = [&](const Vector3& p) -> std::int64_t {
        std::int64_t ix = (std::int64_t)std::floor(p.x * inv);
        std::int64_t iy = (std::int64_t)std::floor(p.y * inv);
        std::int64_t iz = (std::int64_t)std::floor(p.z * inv);
        return ((ix & 0x1FFFFF) << 42) | ((iy & 0x1FFFFF) << 21) | (iz & 0x1FFFFF);
    };
    for (int pass = 0; pass < iterations; ++pass) {
        std::unordered_map<std::int64_t, std::vector<int>> grid;
        grid.reserve(pos.size() * 2);
        for (int i = 0; i < (int)pos.size(); ++i) grid[key(pos[i])].push_back(i);
        for (int i = 0; i < (int)pos.size(); ++i) {
            std::int64_t bx = (std::int64_t)std::floor(pos[i].x * inv);
            std::int64_t by = (std::int64_t)std::floor(pos[i].y * inv);
            std::int64_t bz = (std::int64_t)std::floor(pos[i].z * inv);
            for (std::int64_t dx = -1; dx <= 1; ++dx)
                for (std::int64_t dy = -1; dy <= 1; ++dy)
                    for (std::int64_t dz = -1; dz <= 1; ++dz) {
                        std::int64_t k = (((bx + dx) & 0x1FFFFF) << 42)
                                       | (((by + dy) & 0x1FFFFF) << 21)
                                       | ((bz + dz) & 0x1FFFFF);
                        auto it = grid.find(k);
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            if (j <= i) continue;
                            Vector3 delta = pos[j] - pos[i];
                            real len = delta.magnitude();
                            if (len >= diam) continue;
                            Vector3 n = len > (real)1e-9 ? delta * (((real)1) / len)
                                                         : Vector3(1, 0, 0);
                            real push = (diam - len) * (real)0.5 + (real)1e-9;
                            pos[i] -= n * push;
                            pos[j] += n * push;
                        }
                    }
        }
    }
}

// ============================================================= 3. PBFFluid (PBD)
// Position-Based Fluids: predict positions under gravity, then iteratively project
// them onto the incompressibility constraint C_i = ρ_i/ρ0 − 1 = 0 using the
// per-particle multiplier λ_i and the spiky-kernel constraint gradients.
class PBFFluid {
public:
    real h = (real)0.1;             // smoothing radius
    real mass = 1;                  // per-particle mass
    real rho0 = 1000;               // rest density (measured in build())
    int solverIters = 4;            // constraint-projection iterations per step
    Vector3 gravity = Vector3(0, -9.81, 0);
    real cfmEps = (real)1e2;        // constraint-force-mixing relaxation
    real xsph   = (real)0.0;        // XSPH viscosity coefficient
    // artificial-pressure (tensile-instability) term  −k (W(r)/W(Δq))ⁿ
    real sCorrK = (real)0.0, sCorrDq = (real)0.2, sCorrN = 4;
    Vector3 boundsMin, boundsMax;

    std::vector<Vector3> pos, vel, pred;
    std::vector<real> lambda, rho;

    int nParticles() const { return (int)pos.size(); }

    // --- kernels (3D) ---
    real poly6(real r2) const { if (r2 >= h2) return 0; real u = h2 - r2; return cPoly6 * u * u * u; }
    Vector3 spikyGrad(const Vector3& rij, real r) const {
        if (r >= h || r < (real)1e-9) return Vector3();
        return rij * (cSpiky * (h - r) * (h - r) / r);
    }
    void setKernelConstants() {
        h2 = h * h;
        cPoly6 = (real)315.0 / ((real)64.0 * real_pi * real_pow(h, 9));
        cSpiky = (real)-45.0 / (real_pi * real_pow(h, 6));
    }

    // Fill an axis-aligned block of fluid and measure ρ0 from the packing, so the
    // resting lattice starts at ~zero constraint (as sph.h's build() does).
    void build(const Vector3& cmin, const Vector3& cmax,
               const Vector3& fmin, const Vector3& fmax, real spacing) {
        boundsMin = cmin; boundsMax = cmax; setKernelConstants();
        pos.clear();
        for (real x = fmin.x; x <= fmax.x + (real)1e-9; x += spacing)
            for (real y = fmin.y; y <= fmax.y + (real)1e-9; y += spacing)
                for (real z = fmin.z; z <= fmax.z + (real)1e-9; z += spacing)
                    pos.push_back(Vector3(x, y, z));
        vel.assign(pos.size(), Vector3());
        pred = pos; lambda.assign(pos.size(), 0); rho.assign(pos.size(), rho0);
        buildGrid(pos); computeDensity(pos);
        real sum = 0; for (real d : rho) sum += d;
        if (!rho.empty()) rho0 = sum / rho.size();
    }

    void step(real dt) {
        setKernelConstants();
        for (size_t i = 0; i < pos.size(); ++i) {          // predict
            vel[i] += gravity * dt;
            pred[i] = pos[i] + vel[i] * dt;
            clampToBounds(pred[i]);
        }
        for (int it = 0; it < solverIters; ++it) {
            buildGrid(pred);
            computeDensity(pred);
            for (int i = 0; i < (int)pred.size(); ++i) {   // λ_i
                real C = rho[i] / rho0 - 1;
                Vector3 gradI; real sumSq = 0;
                Vector3 pi = pred[i];
                forNeighbours(pi, [&](int j) {
                    if (j == i) return;
                    Vector3 rij = pi - pred[j]; real r = rij.magnitude();
                    Vector3 g = spikyGrad(rij, r) * (mass / rho0);
                    gradI += g; sumSq += g * g;            // ∇_{p_j}C_i = −g ⇒ |·|²=|g|²
                });
                sumSq += gradI * gradI;
                lambda[i] = -C / (sumSq + cfmEps);
            }
            for (int i = 0; i < (int)pred.size(); ++i) {   // Δp_i
                Vector3 dp, pi = pred[i]; real li = lambda[i];
                forNeighbours(pi, [&](int j) {
                    if (j == i) return;
                    Vector3 rij = pi - pred[j]; real r = rij.magnitude();
                    real sc = 0;
                    if (sCorrK > 0) {
                        real wq = poly6(sCorrDq * sCorrDq * h2);   // W at Δq·h
                        if (wq > 1e-12) sc = -sCorrK * real_pow(poly6(r * r) / wq, sCorrN);
                    }
                    dp += spikyGrad(rij, r) * ((li + lambda[j] + sc) * mass / rho0);
                });
                pred[i] += dp;
                clampToBounds(pred[i]);
            }
        }
        for (size_t i = 0; i < pos.size(); ++i) vel[i] = (pred[i] - pos[i]) * (((real)1) / dt);
        if (xsph > 0) applyXSPH();
        pos = pred;
    }

    // Mean poly6 density over all particles at the current positions.
    real meanDensity() {
        buildGrid(pos); computeDensity(pos);
        real s = 0; for (real d : rho) s += d;
        return rho.empty() ? 0 : s / rho.size();
    }

private:
    real h2 = (real)0.01, cPoly6 = 1, cSpiky = 1;
    int gnx = 1, gny = 1, gnz = 1; Vector3 gorigin; real gcell = (real)0.1;
    std::vector<std::vector<int>> cells;

    void clampToBounds(Vector3& p) const {
        const real e = (real)1e-4;
        for (int a = 0; a < 3; ++a) {
            real lo = boundsMin[a] + e, hi = boundsMax[a] - e;
            if (p[a] < lo) { if (a == 0) p.x = lo; else if (a == 1) p.y = lo; else p.z = lo; }
            if (p[a] > hi) { if (a == 0) p.x = hi; else if (a == 1) p.y = hi; else p.z = hi; }
        }
    }
    int clampi(int v, int hi) const { return v < 0 ? 0 : (v > hi ? hi : v); }
    int cellIndex(const Vector3& p) const {
        int cx = clampi((int)((p.x - gorigin.x) / gcell), gnx - 1);
        int cy = clampi((int)((p.y - gorigin.y) / gcell), gny - 1);
        int cz = clampi((int)((p.z - gorigin.z) / gcell), gnz - 1);
        return (cx * gny + cy) * gnz + cz;
    }
    void buildGrid(const std::vector<Vector3>& P) {
        gcell = h; gorigin = boundsMin;
        gnx = std::max(1, (int)std::ceil((boundsMax.x - boundsMin.x) / gcell) + 1);
        gny = std::max(1, (int)std::ceil((boundsMax.y - boundsMin.y) / gcell) + 1);
        gnz = std::max(1, (int)std::ceil((boundsMax.z - boundsMin.z) / gcell) + 1);
        cells.assign((size_t)gnx * gny * gnz, {});
        for (int i = 0; i < (int)P.size(); ++i) cells[cellIndex(P[i])].push_back(i);
    }
    template <class F> void forNeighbours(const Vector3& p, F&& fn) const {
        int cx = clampi((int)((p.x - gorigin.x) / gcell), gnx - 1);
        int cy = clampi((int)((p.y - gorigin.y) / gcell), gny - 1);
        int cz = clampi((int)((p.z - gorigin.z) / gcell), gnz - 1);
        for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) for (int dz = -1; dz <= 1; ++dz) {
            int x = cx + dx, y = cy + dy, z = cz + dz;
            if (x < 0 || y < 0 || z < 0 || x >= gnx || y >= gny || z >= gnz) continue;
            for (int j : cells[(x * gny + y) * gnz + z]) fn(j);
        }
    }
    void computeDensity(const std::vector<Vector3>& P) {
        rho.assign(P.size(), 0);
        for (int i = 0; i < (int)P.size(); ++i) {
            real d = 0; Vector3 pi = P[i];
            forNeighbours(pi, [&](int j) { real r2 = (pi - P[j]).squareMagnitude(); if (r2 < h2) d += mass * poly6(r2); });
            rho[i] = d;
        }
    }
    void applyXSPH() {
        std::vector<Vector3> corr(pos.size());
        buildGrid(pred); computeDensity(pred);
        for (int i = 0; i < (int)pred.size(); ++i) {
            Vector3 c, pi = pred[i], vi = vel[i];
            forNeighbours(pi, [&](int j) {
                if (j == i) return;
                real r2 = (pi - pred[j]).squareMagnitude();
                if (r2 < h2 && rho[j] > 1e-9) c += (vel[j] - vi) * (mass / rho[j] * poly6(r2));
            });
            corr[i] = c * xsph;
        }
        for (size_t i = 0; i < vel.size(); ++i) vel[i] += corr[i];
    }
};

} // namespace phys
