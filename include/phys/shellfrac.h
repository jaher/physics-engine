// Shell fracture of a thin-walled hollow cylinder — a drum/barrel/pipe blowing apart
// into small curved metal fragments. The thin wall is diced into a jittered grid of
// nCirc × nRows segments (each a small curved hexahedron spanning an arc × a height
// band × the wall thickness), optionally with the two end caps diced into wedges.
// Each piece is a convex ConvexCell with exact mass properties (reusing voronoi.h),
// so they collide and tumble like real shrapnel.
#pragma once
#include "voronoi.h"
#include "glassfrac.h"   // V2 + gClipHalf (2-D half-plane clipping)
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

namespace phys {

// build a hexahedron (8 corners; 0-3 one quad, 4-7 the opposite quad, i↔i+4 paired)
// into a ConvexCell, orienting faces outward robustly before the mass-property solve.
inline bool buildHex(const std::array<Vector3, 8>& p, ConvexCell& cell) {
    cell = ConvexCell();
    for (auto& v : p) cell.verts.push_back(v);
    const int F[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 6, 5}};
    for (auto& f : F) { cell.tris.push_back({f[0], f[1], f[2]}); cell.tris.push_back({f[0], f[2], f[3]}); }
    Vector3 mc(0, 0, 0); for (auto& v : cell.verts) mc += v; mc = mc * (1.0 / 8.0);
    cell.triN.assign(cell.tris.size(), Vector3());
    for (size_t ti = 0; ti < cell.tris.size(); ti++) {
        auto& tr = cell.tris[ti];
        Vector3 a = cell.verts[tr[0]], b = cell.verts[tr[1]], c = cell.verts[tr[2]];
        Vector3 gn = (b - a) % (c - a); Vector3 tc = (a + b + c) * (1.0 / 3.0);
        if (gn * (tc - mc) < 0) { std::swap(tr[1], tr[2]); gn = gn * -1; }
        real l = gn.magnitude(); cell.triN[ti] = l > 1e-12 ? gn * (1.0 / l) : Vector3(0, 1, 0);
    }
    cellMassProps(cell);
    if (cell.volume < 1e-11) return false;
    cell.site = cell.com;
    for (auto& v : cell.verts) v -= cell.com;
    return true;
}

// Fracture a hollow cylinder (mean radius R, half-height Hh, wall thickness wt, axis
// along y) into small shell pieces: nCirc × nRows wall segments plus, if `caps`, the
// top & bottom end plates diced into nCirc wedges each. Boundaries are jittered.
inline std::vector<ConvexCell> cylinderShellFracture(real R, real Hh, real wt, int nCirc, int nRows,
                                                     bool caps = true, unsigned seed = 7u) {
    unsigned s = seed ? seed : 1u; auto rr = [&]() { s = s * 1103515245u + 12345u; return (real)(((s >> 16) & 0x7fff) / 32767.0); };
    const real ri = R - wt * 0.5, ro = R + wt * 0.5;
    std::vector<real> ang(nCirc + 1), hgt(nRows + 1);
    for (int j = 0; j <= nCirc; j++) ang[j] = 2 * M_PI * j / nCirc;
    for (int i = 0; i <= nRows; i++) hgt[i] = -Hh + 2 * Hh * i / nRows;
    for (int j = 1; j < nCirc; j++) ang[j] += (rr() - 0.5) * (2 * M_PI / nCirc) * 0.45;   // jitter interior seams
    for (int i = 1; i < nRows; i++) hgt[i] += (rr() - 0.5) * (2 * Hh / nRows) * 0.45;
    auto P = [&](real a, real y, real r) { return Vector3(r * std::cos(a), y, r * std::sin(a)); };
    std::vector<ConvexCell> cells;
    for (int i = 0; i < nRows; i++) for (int j = 0; j < nCirc; j++) {           // wall segments
        real a0 = ang[j], a1 = ang[j + 1], y0 = hgt[i], y1 = hgt[i + 1];
        std::array<Vector3, 8> p = {P(a0, y0, ri), P(a1, y0, ri), P(a1, y1, ri), P(a0, y1, ri),
                                    P(a0, y0, ro), P(a1, y0, ro), P(a1, y1, ro), P(a0, y1, ro)};
        ConvexCell c; if (buildHex(p, c)) cells.push_back(std::move(c));
    }
    if (caps) for (int cap = 0; cap < 2; cap++) {                                // end plates → wedges
        real yc = cap ? Hh : -Hh, y0 = yc - wt * 0.5, y1 = yc + wt * 0.5, rh = R * 0.14;
        for (int j = 0; j < nCirc; j++) { real a0 = ang[j], a1 = ang[j + 1];
            std::array<Vector3, 8> p = {P(a0, y0, rh), P(a1, y0, rh), P(a1, y0, ro), P(a0, y0, ro),
                                        P(a0, y1, rh), P(a1, y1, rh), P(a1, y1, ro), P(a0, y1, ro)};
            ConvexCell c; if (buildHex(p, c)) cells.push_back(std::move(c)); }
    }
    return cells;
}

// map a 2-D wall polygon (in arc-length s × height v) onto the cylinder and extrude it
// through the wall as a shell shard (a ConvexCell, recentred about its COM). The polygon
// is placed on the tangent plane at the cell centre (each shard spans a small arc, so
// the flat approximation is negligible) — planar faces keep the mass solve robust.
inline bool buildShellPiece(const std::vector<V2>& poly, real R, real ri, real ro, ConvexCell& cell) {
    int N = (int)poly.size(); if (N < 3) return false;
    real sc = 0, vc = 0; for (auto& p : poly) { sc += p.u; vc += p.v; } sc /= N; vc /= N;
    real uc = sc / R; Vector3 Pc(R * std::cos(uc), vc, R * std::sin(uc));
    Vector3 tang(-std::sin(uc), 0, std::cos(uc)), up(0, 1, 0), nrm(std::cos(uc), 0, std::sin(uc));   // local frame on the wall
    real hw = (ro - ri) * 0.5;
    cell = ConvexCell();
    for (auto& p : poly) { Vector3 b = Pc + tang * (p.u - sc) + up * (p.v - vc); cell.verts.push_back(b + nrm * hw); }   // outer ring 0..N-1
    for (auto& p : poly) { Vector3 b = Pc + tang * (p.u - sc) + up * (p.v - vc); cell.verts.push_back(b - nrm * hw); }   // inner ring N..2N-1
    for (int i = 1; i < N - 1; i++) cell.tris.push_back({0, i, i + 1});                          // outer face
    for (int i = 1; i < N - 1; i++) cell.tris.push_back({N, N + i, N + i + 1});                  // inner face
    for (int i = 0; i < N; i++) { int j = (i + 1) % N; cell.tris.push_back({i, j, N + j}); cell.tris.push_back({i, N + j, N + i}); }   // sides
    Vector3 mc(0, 0, 0); for (auto& v : cell.verts) mc += v; mc = mc * (1.0 / cell.verts.size());
    cell.triN.assign(cell.tris.size(), Vector3());
    for (size_t ti = 0; ti < cell.tris.size(); ti++) {
        auto& tr = cell.tris[ti];
        Vector3 a = cell.verts[tr[0]], b = cell.verts[tr[1]], c = cell.verts[tr[2]];
        Vector3 gn = (b - a) % (c - a); Vector3 tc = (a + b + c) * (1.0 / 3.0);
        if (gn * (tc - mc) < 0) { std::swap(tr[1], tr[2]); gn = gn * -1; }
        real l = gn.magnitude(); cell.triN[ti] = l > 1e-12 ? gn * (1.0 / l) : Vector3(0, 1, 0);
    }
    cellMassProps(cell);
    if (cell.volume < 1e-11) return false;
    cell.site = cell.com;
    for (auto& v : cell.verts) v -= cell.com;
    return true;
}

// Voronoi-style shell fracture: shatter a hollow cylinder wall into IRREGULAR convex
// metal shards. Seeds are a jittered nCirc×nRows grid on the unrolled wall (arc-length
// s × height v) — jittered for irregular shapes but size-bounded so every shard spans a
// small arc; a 2-D Voronoi is computed there with s PERIODIC (seeds replicated at s ± 2πR
// so cells tile seamlessly across the seam), and each cell is mapped back onto the
// cylinder and extruded through the wall thickness.
inline std::vector<ConvexCell> cylinderShellVoronoi(real R, real Hh, real wt, int nCirc, int nRows, unsigned seed = 7u) {
    unsigned s = seed ? seed : 1u; auto rr = [&]() { s = s * 1103515245u + 12345u; return (real)(((s >> 16) & 0x7fff) / 32767.0); };
    const real circ = 2 * M_PI * R, ri = R - wt * 0.5, ro = R + wt * 0.5, hc = circ * 0.5;
    std::vector<V2> seeds;
    for (int cj = 0; cj < nCirc; cj++) for (int ci = 0; ci < nRows; ci++)
        seeds.push_back({(cj + 0.5 + (rr() - 0.5) * 0.8) / nCirc * circ, -Hh + (ci + 0.5 + (rr() - 0.5) * 0.8) / nRows * 2 * Hh});
    int nSeeds = (int)seeds.size();
    std::vector<ConvexCell> cells;
    for (int i = 0; i < nSeeds; i++) {
        V2 si = seeds[i];
        std::vector<V2> poly = {{si.u - hc, -Hh}, {si.u + hc, -Hh}, {si.u + hc, Hh}, {si.u - hc, Hh}};   // ≤ half-circumference cell
        for (int j = 0; j < nSeeds; j++) for (int rep = -1; rep <= 1; rep++) {
            if (j == i && rep == 0) continue;
            real su = seeds[j].u + rep * circ, sv = seeds[j].v;
            if (std::fabs(su - si.u) > hc + 1e-3) continue;
            real nu = su - si.u, nv = sv - si.v, d = 0.5 * (nu * (si.u + su) + nv * (si.v + sv));         // perpendicular bisector, keep si side
            poly = gClipHalf(poly, nu, nv, d);
        }
        poly = gClipHalf(poly, 0, 1, Hh); poly = gClipHalf(poly, 0, -1, Hh);
        if (poly.size() < 3) continue;
        real area2D = 0; for (int k = 0; k < (int)poly.size(); k++) { auto& p0 = poly[k]; auto& p1 = poly[(k + 1) % poly.size()]; area2D += p0.u * p1.v - p1.u * p0.v; }
        area2D = std::fabs(area2D) * 0.5; if (area2D < 1e-4) continue;                                // drop degenerate slivers
        ConvexCell cell; if (!buildShellPiece(poly, R, ri, ro, cell)) continue;
        real ev = area2D * wt;                                                                        // a shell piece's volume ≈ (unrolled area)·thickness
        if (cell.volume < 0.6 * ev || cell.volume > 1.5 * ev) continue;                              // reject any piece whose mass-solve went bad
        if (cell.inertiaUnit.data[0] <= 0 || cell.inertiaUnit.data[4] <= 0 || cell.inertiaUnit.data[8] <= 0) continue;
        cells.push_back(std::move(cell));
    }
    return cells;
}

} // namespace phys
