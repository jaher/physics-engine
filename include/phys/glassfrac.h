// Impact ("bullet-hole") fracture of a thin glass pane. Real glass struck by a
// projectile breaks into the classic spider-web: a small punched-out hole, RADIAL
// cracks spoking outward, and CONCENTRIC rings joining them — fine shards packed
// around the entry point, coarser plates toward the rim. We build exactly that: a
// jittered polar grid of cracks centred on the impact, each cell extruded to the
// glass thickness into a thin convex-polyhedron fragment (reusing ConvexCell + the
// exact mass-property solver from voronoi.h). Cells are clipped to the pane rect.
#pragma once
#include "voronoi.h"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

namespace phys {

struct V2 { real u, v; };

// clip a convex polygon by the half-plane  n·x ≤ d  (2-D)
inline std::vector<V2> gClipHalf(const std::vector<V2>& poly, real nu, real nv, real d) {
    std::vector<V2> out; int N = (int)poly.size();
    for (int i = 0; i < N; i++) {
        V2 A = poly[i], B = poly[(i + 1) % N];
        real da = nu * A.u + nv * A.v - d, db = nu * B.u + nv * B.v - d;
        if (da <= 0) out.push_back(A);
        if ((da < 0) != (db < 0)) { real t = da / (da - db); out.push_back({A.u + t * (B.u - A.u), A.v + t * (B.v - A.v)}); }
    }
    return out;
}
inline std::vector<V2> gClipRect(std::vector<V2> p, real w, real h) {
    p = gClipHalf(p, 1, 0, w); p = gClipHalf(p, -1, 0, w); p = gClipHalf(p, 0, 1, h); p = gClipHalf(p, 0, -1, h); return p;
}

// extrude a 2-D convex polygon by ±t/2 along local z into a ConvexCell (recentred
// about its COM). Faces are oriented outward robustly (convex shape) before the
// mass-property solve, so winding never matters to the caller.
inline bool gBuildPrism(const std::vector<V2>& poly, real t, ConvexCell& cell) {
    int N = (int)poly.size(); if (N < 3) return false;
    cell = ConvexCell();
    for (int i = 0; i < N; i++) cell.verts.push_back(Vector3(poly[i].u, poly[i].v, t * 0.5));   // front 0..N-1
    for (int i = 0; i < N; i++) cell.verts.push_back(Vector3(poly[i].u, poly[i].v, -t * 0.5));  // back  N..2N-1
    for (int i = 1; i < N - 1; i++) cell.tris.push_back({0, i, i + 1});                          // front fan
    for (int i = 1; i < N - 1; i++) cell.tris.push_back({N, N + i, N + i + 1});                  // back fan
    for (int i = 0; i < N; i++) { int j = (i + 1) % N; cell.tris.push_back({i, j, N + j}); cell.tris.push_back({i, N + j, N + i}); }
    Vector3 mc(0, 0, 0); for (auto& v : cell.verts) mc += v; mc = mc * (1.0 / cell.verts.size());
    cell.triN.assign(cell.tris.size(), Vector3());
    for (size_t ti = 0; ti < cell.tris.size(); ti++) {
        auto& tr = cell.tris[ti];
        Vector3 a = cell.verts[tr[0]], b = cell.verts[tr[1]], c = cell.verts[tr[2]];
        Vector3 gn = (b - a) % (c - a); Vector3 tc = (a + b + c) * (1.0 / 3.0);
        if (gn * (tc - mc) < 0) { std::swap(tr[1], tr[2]); gn = gn * -1; }                        // orient outward
        real l = gn.magnitude(); cell.triN[ti] = l > 1e-12 ? gn * (1.0 / l) : Vector3(0, 0, 1);
    }
    cellMassProps(cell);
    if (cell.volume < 1e-9) return false;
    cell.site = cell.com;
    for (auto& v : cell.verts) v -= cell.com;
    return true;
}

// Fracture a pane (half-widths w,h; thickness t) about the impact (iu,iv). `nSpokes`
// radial cracks × `nRings` concentric rings (radii grow geometrically → fine at the
// hole, coarse at the rim); a hole of radius `holeR` is punched out. Spoke angles and
// ring radii are jittered so the web looks natural.
inline std::vector<ConvexCell> glassFracture(real w, real h, real t, real iu, real iv,
        int nSpokes = 24, int nRings = 8, real holeR = 0.05, unsigned seed = 7u) {
    unsigned s = seed ? seed : 1u; auto rr = [&]() { s = s * 1103515245u + 12345u; return (real)(((s >> 16) & 0x7fff) / 32767.0); };
    real maxR = 0;
    for (real cu : {-w, w}) for (real cv : {-h, h}) maxR = std::max(maxR, (real)std::hypot((double)(cu - iu), (double)(cv - iv)));
    maxR *= 1.06;
    std::vector<real> r(nRings + 1); real growth = std::pow(maxR / holeR, 1.0 / nRings);
    for (int k = 0; k <= nRings; k++) r[k] = holeR * std::pow(growth, (real)k);
    std::vector<real> th(nSpokes + 1);
    for (int j = 0; j < nSpokes; j++) th[j] = 2 * M_PI * j / nSpokes + (rr() - 0.5) * (2 * M_PI / nSpokes) * 0.55;
    std::sort(th.begin(), th.begin() + nSpokes); th[nSpokes] = th[0] + 2 * M_PI;
    // shared, jittered ring radius per (ring, spoke) so adjacent cells tile exactly
    std::vector<std::vector<real>> Rg(nRings + 1, std::vector<real>(nSpokes + 1));
    for (int k = 0; k <= nRings; k++) { for (int j = 0; j < nSpokes; j++) Rg[k][j] = r[k] * (k == 0 ? 1.0 : (1.0 + (rr() - 0.5) * 0.28)); Rg[k][nSpokes] = Rg[k][0]; }
    auto pt = [&](int k, int j) { return V2{iu + Rg[k][j] * std::cos(th[j]), iv + Rg[k][j] * std::sin(th[j])}; };
    std::vector<ConvexCell> cells;
    for (int k = 0; k < nRings; k++) for (int j = 0; j < nSpokes; j++) {
        std::vector<V2> quad = {pt(k, j), pt(k + 1, j), pt(k + 1, j + 1), pt(k, j + 1)};
        std::vector<V2> poly = gClipRect(quad, w, h);
        ConvexCell cell; if (gBuildPrism(poly, t, cell)) cells.push_back(std::move(cell));
    }
    return cells;
}

} // namespace phys
