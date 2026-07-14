// Voronoi fracture. Scatter seed points through a box and each fragment becomes
// the Voronoi cell of its seed — the set of points closer to that seed than any
// other, clipped to the box. Every cell is a convex polyhedron with irregular
// angular faces (not a parallelepiped), so a shattered solid looks like real
// cracked stone. Cells are produced as a vertex set (all GJK/EPA collision needs),
// triangulated faces (for rendering), and exact mass properties (volume, centre of
// mass, and inertia tensor of the polyhedron at unit density).
#pragma once
#include "core.h"
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

namespace phys {

struct ConvexCell {
    Vector3 site;                                 // the seed point
    std::vector<Vector3> verts;                   // hull vertices (local space, about the COM)
    std::vector<std::array<int, 3>> tris;         // triangulated faces → indices into verts
    std::vector<Vector3> triN;                    // per-triangle outward normal (flat shading)
    Vector3 com; real volume = 0;                 // mass properties (about the COM) at unit density
    Matrix3 inertiaUnit;
};

// exact mass properties of a closed triangle mesh (outward CCW), unit density.
// Uses the tetra-from-origin covariance sum, then shifts to the centre of mass.
inline void cellMassProps(ConvexCell& c) {
    real vol = 0; Vector3 fm(0, 0, 0);                     // volume + first moment ∫r dV
    // covariance C = ∫ r rᵀ dV about the origin
    real C[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (auto& t : c.tris) {
        Vector3 a = c.verts[t[0]], b = c.verts[t[1]], d = c.verts[t[2]];
        real det = a * (b % d);                            // 6·(signed tet volume)
        vol += det;
        fm += (a + b + d) * det;                           // ×1/24 below
        // Ĉ_ref = (1/120)[[2,1,1],[1,2,1],[1,1,2]];  C_tet = det · M Ĉ Mᵀ,  M=[a b d]
        Vector3 col[3] = {a, b, d};
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
            real ai[3] = {col[0][i], col[1][i], col[2][i]}, aj[3] = {col[0][j], col[1][j], col[2][j]};
            real s = 0;
            for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) s += (p == q ? 2.0 : 1.0) * ai[p] * aj[q];
            C[i][j] += det * s / 120.0;
        }
    }
    real V = vol / 6.0; c.volume = V;
    Vector3 com = (V > 1e-12) ? fm * (1.0 / (24.0 * V)) : Vector3();
    c.com = com;
    // inertia about origin: I = tr(C)·E − C ; then parallel-axis shift origin→com
    real tr = C[0][0] + C[1][1] + C[2][2];
    Matrix3 Io; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Io.data[i * 3 + j] = (i == j ? tr : 0) - C[i][j];
    real d2 = com * com; Matrix3 shift;
    real dd[3] = {com.x, com.y, com.z};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) shift.data[i * 3 + j] = V * ((i == j ? d2 : 0) - dd[i] * dd[j]);
    for (int k = 0; k < 9; k++) c.inertiaUnit.data[k] = Io.data[k] - shift.data[k];
}

// intersection point of three planes nᵢ·x = dᵢ ; returns false if near-parallel.
inline bool planeTriple(const Vector3& n1, real d1, const Vector3& n2, real d2,
                        const Vector3& n3, real d3, Vector3& out) {
    Vector3 c23 = n2 % n3; real den = n1 * c23;
    if (std::fabs(den) < 1e-9) return false;
    out = (c23 * d1 + (n3 % n1) * d2 + (n1 % n2) * d3) * (1.0 / den);
    return true;
}

// build the Voronoi cells of `sites` inside the box [-half, half]. Each cell is
// bounded by the 6 box planes plus the perpendicular bisectors to its nearest
// `kNear` neighbours. Cells whose interior is degenerate are skipped.
inline std::vector<ConvexCell> voronoiFracture(const Vector3& half, const std::vector<Vector3>& sites, int kNear = 18) {
    struct Pl { Vector3 n; real d; };
    std::vector<ConvexCell> cells;
    const real eps = 1e-5 * (half.x + half.y + half.z);
    for (size_t i = 0; i < sites.size(); i++) {
        Vector3 s = sites[i];
        std::vector<Pl> P = {{Vector3(1, 0, 0), half.x}, {Vector3(-1, 0, 0), half.x},
                             {Vector3(0, 1, 0), half.y}, {Vector3(0, -1, 0), half.y},
                             {Vector3(0, 0, 1), half.z}, {Vector3(0, 0, -1), half.z}};
        // nearest neighbours → bisector planes (n·x ≤ d keeps the seed side)
        std::vector<std::pair<real, int>> nb;
        for (size_t j = 0; j < sites.size(); j++) if (j != i) nb.push_back({(sites[j] - s).squareMagnitude(), (int)j});
        std::sort(nb.begin(), nb.end());
        int kk = std::min((int)nb.size(), kNear);
        for (int m = 0; m < kk; m++) { Vector3 dir = sites[nb[m].second] - s; real len = dir.magnitude(); if (len < 1e-9) continue;
            Vector3 n = dir * (1.0 / len); Vector3 mid = (s + sites[nb[m].second]) * 0.5;
            P.push_back({n, n * mid});                     // bisector half-space, world frame (box-centred)
        }
        // enumerate feasible triple-plane vertices (local frame, origin at the seed)
        std::vector<Vector3> V; int np = (int)P.size();
        for (int a = 0; a < np; a++) for (int b = a + 1; b < np; b++) for (int cc = b + 1; cc < np; cc++) {
            Vector3 x; if (!planeTriple(P[a].n, P[a].d, P[b].n, P[b].d, P[cc].n, P[cc].d, x)) continue;
            bool ok = true; for (int k = 0; k < np; k++) if (P[k].n * x > P[k].d + eps) { ok = false; break; }
            if (!ok) continue;
            bool dup = false; for (auto& v : V) if ((v - x).squareMagnitude() < eps * eps) { dup = true; break; }
            if (!dup) V.push_back(x);
        }
        if (V.size() < 4) continue;                        // degenerate cell
        ConvexCell cell; cell.site = s; cell.verts = V;
        // faces: verts lying on each plane, ordered by angle, fan-triangulated
        for (auto& pl : P) {
            std::vector<int> on;
            for (int k = 0; k < (int)V.size(); k++) if (std::fabs(pl.n * V[k] - pl.d) < eps) on.push_back(k);
            if ((int)on.size() < 3) continue;
            Vector3 cen(0, 0, 0); for (int k : on) cen += V[k]; cen = cen * (1.0 / on.size());
            Vector3 u = (V[on[0]] - cen); real ul = u.magnitude(); if (ul < 1e-9) continue; u = u * (1.0 / ul);
            Vector3 w = pl.n % u;
            std::sort(on.begin(), on.end(), [&](int A, int B) {
                Vector3 pa = V[A] - cen, pb = V[B] - cen;
                return std::atan2(pa * w, pa * u) < std::atan2(pb * w, pb * u); });
            for (int k = 1; k + 1 < (int)on.size(); k++) {
                std::array<int, 3> tri = {on[0], on[k], on[k + 1]};
                Vector3 fn = (V[tri[1]] - V[tri[0]]) % (V[tri[2]] - V[tri[0]]);
                if (fn * pl.n < 0) std::swap(tri[1], tri[2]);   // outward winding
                cell.tris.push_back(tri); cell.triN.push_back(pl.n);
            }
        }
        if (cell.tris.size() < 4) continue;
        cellMassProps(cell);
        if (cell.volume < 1e-9) continue;
        // vertices are in absolute (box-centred) coordinates; recentre them about the
        // cell's COM so the body sits at the COM and the hull is COM-local.
        cell.site = cell.com;                              // cell centroid in the box frame
        for (auto& v : cell.verts) v -= cell.com;
        cells.push_back(std::move(cell));
    }
    return cells;
}

} // namespace phys
