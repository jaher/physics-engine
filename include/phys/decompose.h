// Approximate convex decomposition (V-HACD-lite) for a 2D polygon. Physics engines
// only collide *convex* shapes cheaply, so a concave outline (an L, a star, a
// letter) has to be broken into a few convex pieces first. We do it exactly rather
// than by voxel approximation: find a reflex (concave) vertex, cut the polygon with
// an interior diagonal from it, and recurse on the two halves until every piece is
// convex. The diagonal test is O'Rourke's InCone + no-edge-crossing check, so the
// pieces tile the input exactly (their areas sum to the original). An L-shape comes
// out as two convex quads; each extra reflex feature adds roughly one more cut.
// Polygons live in the XY plane as `Vector3` (z ignored) to compose with the engine.
#pragma once
#include "core.h"
#include <vector>
#include <algorithm>

namespace phys {

using Polygon2 = std::vector<Vector3>;

// twice the signed area of triangle (a,b,c) in XY; > 0 ⇔ c is left of a→b (CCW turn)
inline real area2(const Vector3& a, const Vector3& b, const Vector3& c) {
    return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}
inline bool left  (const Vector3& a, const Vector3& b, const Vector3& c) { return area2(a, b, c) >  0; }
inline bool leftOn(const Vector3& a, const Vector3& b, const Vector3& c) { return area2(a, b, c) >= 0; }
inline bool collinear(const Vector3& a, const Vector3& b, const Vector3& c) { return area2(a, b, c) == 0; }

// twice the signed area of a polygon (shoelace); > 0 ⇔ counter-clockwise
inline real polygonArea2(const Polygon2& p) {
    real s = 0; int n = (int)p.size();
    for (int i = 0; i < n; i++) { const Vector3& a = p[i]; const Vector3& b = p[(i + 1) % n]; s += a.x * b.y - b.x * a.y; }
    return s;
}
inline real polygonArea(const Polygon2& p) { return real_abs(polygonArea2(p)) * (real)0.5; }

// Is vertex i a reflex (concave) corner of a CCW polygon? Interior angle > 180°
// shows up as a clockwise turn at i.
inline bool isReflex(const Polygon2& p, int i) {
    int n = (int)p.size();
    return area2(p[(i - 1 + n) % n], p[i], p[(i + 1) % n]) < 0;
}

// A polygon is convex if no vertex is reflex (allowing collinear runs).
inline bool isConvex(const Polygon2& p) {
    int n = (int)p.size();
    if (n < 3) return false;
    for (int i = 0; i < n; i++) if (isReflex(p, i)) return false;
    return true;
}

// --- segment intersection (O'Rourke, "Computational Geometry in C") -------------
inline bool xorb(bool x, bool y) { return x != y; }
inline bool intersectProp(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    if (collinear(a, b, c) || collinear(a, b, d) || collinear(c, d, a) || collinear(c, d, b)) return false;
    return xorb(left(a, b, c), left(a, b, d)) && xorb(left(c, d, a), left(c, d, b));
}
inline bool between(const Vector3& a, const Vector3& b, const Vector3& c) { // c collinear-with & between a,b
    if (!collinear(a, b, c)) return false;
    if (a.x != b.x) return (a.x <= c.x && c.x <= b.x) || (a.x >= c.x && c.x >= b.x);
    return (a.y <= c.y && c.y <= b.y) || (a.y >= c.y && c.y >= b.y);
}
inline bool segIntersect(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    if (intersectProp(a, b, c, d)) return true;
    return between(a, b, c) || between(a, b, d) || between(c, d, a) || between(c, d, b);
}

// Does the segment i→j lie in the interior cone at vertex i? (O'Rourke InCone)
inline bool inCone(const Polygon2& p, int i, int j) {
    int n = (int)p.size();
    const Vector3& a = p[i]; const Vector3& a1 = p[(i + 1) % n]; const Vector3& a0 = p[(i - 1 + n) % n]; const Vector3& b = p[j];
    if (leftOn(a, a1, a0))                                  // i is a convex vertex
        return left(a, b, a0) && left(b, a, a1);
    return !(leftOn(a, b, a1) && leftOn(b, a, a0));         // i is reflex
}
// Is i→j a valid internal diagonal: in both cones and crossing no edge?
inline bool diagonal(const Polygon2& p, int i, int j) {
    int n = (int)p.size();
    if (i == j || (i + 1) % n == j || (j + 1) % n == i) return false;   // same / adjacent
    if (!inCone(p, i, j) || !inCone(p, j, i)) return false;
    const Vector3& a = p[i]; const Vector3& b = p[j];
    for (int k = 0; k < n; k++) {
        int k1 = (k + 1) % n;
        if (k == i || k1 == i || k == j || k1 == j) continue;
        if (segIntersect(a, b, p[k], p[k1])) return false;
    }
    return true;
}

// vertices from i to j inclusive, walking forward (cyclically)
inline Polygon2 slicePolygon(const Polygon2& p, int i, int j) {
    int n = (int)p.size(); Polygon2 out;
    for (int k = i;; k = (k + 1) % n) { out.push_back(p[k]); if (k == j) break; }
    return out;
}

// Recurse: cut at the first reflex vertex along an interior diagonal, splitting the
// polygon into two smaller ones, until no reflex vertices remain (convex).
inline void convexDecomposeRec(const Polygon2& poly, std::vector<Polygon2>& out, int depth) {
    int n = (int)poly.size();
    if (n < 3) return;
    if (n == 3 || isConvex(poly) || depth <= 0) { out.push_back(poly); return; }

    int r = -1;
    for (int i = 0; i < n; i++) if (isReflex(poly, i)) { r = i; break; }
    if (r < 0) { out.push_back(poly); return; }            // convex after all

    // pick the best diagonal out of the reflex vertex: prefer cutting to another
    // reflex vertex (removes two concavities at once), else the shortest cut.
    int best = -1; real bestScore = -REAL_MAX;
    for (int j = 0; j < n; j++) {
        if (!diagonal(poly, r, j)) continue;
        Vector3 dv = poly[j] - poly[r];
        real score = (isReflex(poly, j) ? (real)1e12 : (real)0) - dv.squareMagnitude();
        if (score > bestScore) { bestScore = score; best = j; }
    }
    // Fallback (degenerate geometry): any diagonal anywhere keeps recursion progressing.
    if (best < 0) {
        for (int i = 0; i < n && best < 0; i++)
            for (int j = 0; j < n; j++)
                if (diagonal(poly, i, j)) { r = i; best = j; break; }
    }
    if (best < 0) { out.push_back(poly); return; }         // no diagonal — give up cleanly

    convexDecomposeRec(slicePolygon(poly, r, best), out, depth - 1);
    convexDecomposeRec(slicePolygon(poly, best, r), out, depth - 1);
}

// Public entry point: normalise to CCW, then decompose. Returns a set of convex
// pieces that tile `poly` (union area == area of `poly`).
inline std::vector<Polygon2> convexDecompose(Polygon2 poly, int maxDepth = 64) {
    std::vector<Polygon2> out;
    if (poly.size() < 3) return out;
    if (polygonArea2(poly) < 0) std::reverse(poly.begin(), poly.end());   // make CCW
    convexDecomposeRec(poly, out, maxDepth);
    return out;
}

} // namespace phys
