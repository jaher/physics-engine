// Shell fracture of a thin-walled hollow cylinder — a drum/barrel/pipe blowing apart
// into small curved metal fragments. The thin wall is diced into a jittered grid of
// nCirc × nRows segments (each a small curved hexahedron spanning an arc × a height
// band × the wall thickness), optionally with the two end caps diced into wedges.
// Each piece is a convex ConvexCell with exact mass properties (reusing voronoi.h),
// so they collide and tumble like real shrapnel.
#pragma once
#include "voronoi.h"
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

} // namespace phys
