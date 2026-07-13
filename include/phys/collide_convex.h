// Convex & mesh geometry — the shape family every mainstream engine ships and
// this one lacked. Adds, on top of collide_fine.h's sphere/box/plane world:
//   * CollisionConvex   — arbitrary convex hull (point set) on a RigidBody, with
//                         GJK boolean + closest-distance and EPA penetration,
//                         driving a convex–convex contact detector.
//   * CollisionCylinder /
//     CollisionCone     — analytic primitives vs plane / sphere / box.
//   * CollisionTriMesh  — static triangle-soup collider (per-triangle AABB +
//                         optional uniform-grid midphase) vs sphere / box.
// Everything writes standard `Contact`s into a `CollisionData`, same contact
// convention as collide_fine.h (normal points from body[1] toward body[0], the
// direction body[0] is pushed to separate).
#pragma once
#include "collide_fine.h"
#include <vector>
#include <unordered_map>
#include <utility>
#include <cmath>

namespace phys {

// ------------------------------------------------------------- convex primitive
// A convex hull given by local-space vertices, attached to a RigidBody like any
// other CollisionPrimitive. The support() function (farthest vertex along a
// world-space direction) is all GJK/EPA need.
class CollisionConvex : public CollisionPrimitive {
public:
    std::vector<Vector3> vertices;              // local-space point set (the hull)

    // convenience: fill with the 8 corners of a box of the given half-size
    void setBox(const Vector3& h) {
        vertices = {
            Vector3( h.x,  h.y,  h.z), Vector3(-h.x,  h.y,  h.z),
            Vector3( h.x, -h.y,  h.z), Vector3(-h.x, -h.y,  h.z),
            Vector3( h.x,  h.y, -h.z), Vector3(-h.x,  h.y, -h.z),
            Vector3( h.x, -h.y, -h.z), Vector3(-h.x, -h.y, -h.z)};
    }
    // farthest point of the hull along a world-space direction
    Vector3 support(const Vector3& worldDir) const {
        if (vertices.empty()) return getAxis(3);
        Vector3 ld = getTransform().transformInverseDirection(worldDir);   // dir → local
        real best = -REAL_MAX; size_t bi = 0;
        for (size_t i = 0; i < vertices.size(); i++) {
            real d = vertices[i] * ld;
            if (d > best) { best = d; bi = i; }
        }
        return getTransform().transform(vertices[bi]);
    }
};

// ------------------------------------------------------------- round primitives
// Cylinder along its local Y axis: radius r, extent ±halfHeight.
class CollisionCylinder : public CollisionPrimitive {
public:
    real radius = (real)0.5, halfHeight = (real)0.5;
    void endpoints(Vector3& top, Vector3& bot) const {
        Vector3 c = getAxis(3), ax = getAxis(1);
        top = c + ax * halfHeight; bot = c - ax * halfHeight;
    }
    // exact support: cap chosen by axial sign, plus the radial rim offset
    Vector3 support(const Vector3& dir) const {
        Vector3 c = getAxis(3), ax = getAxis(1);
        Vector3 p = c + ax * (dir * ax >= 0 ? halfHeight : -halfHeight);
        Vector3 radial = dir - ax * (dir * ax); real rl = radial.magnitude();
        if (rl > (real)1e-9) p += radial * (radius / rl);
        return p;
    }
};

// Cone along local Y: apex at +Y*halfHeight, base circle (radius r) at -Y*halfHeight.
class CollisionCone : public CollisionPrimitive {
public:
    real radius = (real)0.5, halfHeight = (real)0.5;
    Vector3 apex() const { return getAxis(3) + getAxis(1) * halfHeight; }
    Vector3 baseCentre() const { return getAxis(3) - getAxis(1) * halfHeight; }
    // support = whichever is farther along dir: the apex or the base rim extreme
    Vector3 support(const Vector3& dir) const {
        Vector3 ax = getAxis(1);
        Vector3 ap = apex();
        Vector3 radial = dir - ax * (dir * ax); real rl = radial.magnitude();
        Vector3 base = baseCentre() + (rl > (real)1e-9 ? radial * (radius / rl) : Vector3());
        return (dir * ap >= dir * base) ? ap : base;
    }
};

// ------------------------------------------------------------- GJK / EPA guts
namespace cvx {

// A point on the Minkowski difference A⊖B, keeping the two originating support
// points so witnesses (contact points) can be recovered by barycentric blend.
struct Support { Vector3 v, a, b; };

inline Support mdSupport(const CollisionConvex& A, const CollisionConvex& B, const Vector3& dir) {
    Support s; s.a = A.support(dir); s.b = B.support(-dir); s.v = s.a - s.b; return s;
}

inline Vector3 closestOnSeg(const Vector3& a, const Vector3& b, const Vector3& p) {
    Vector3 ab = b - a; real d = ab * ab;
    real t = d > (real)1e-18 ? ((p - a) * ab) / d : 0;
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    return a + ab * t;
}

// Closest point on triangle abc to p; also returns barycentric (u,v,w) (Ericson).
inline Vector3 closestPtTriangle(const Vector3& p, const Vector3& a, const Vector3& b,
                                 const Vector3& c, real& u, real& v, real& w) {
    Vector3 ab = b - a, ac = c - a, ap = p - a;
    real d1 = ab * ap, d2 = ac * ap;
    if (d1 <= 0 && d2 <= 0) { u = 1; v = 0; w = 0; return a; }
    Vector3 bp = p - b; real d3 = ab * bp, d4 = ac * bp;
    if (d3 >= 0 && d4 <= d3) { u = 0; v = 1; w = 0; return b; }
    real vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) { real t = d1 / (d1 - d3); u = 1 - t; v = t; w = 0; return a + ab * t; }
    Vector3 cp = p - c; real d5 = ab * cp, d6 = ac * cp;
    if (d6 >= 0 && d5 <= d6) { u = 0; v = 0; w = 1; return c; }
    real vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) { real t = d2 / (d2 - d6); u = 1 - t; v = 0; w = t; return a + ac * t; }
    real va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        real t = (d4 - d3) / ((d4 - d3) + (d5 - d6)); u = 0; v = 1 - t; w = t; return b + (c - b) * t;
    }
    real denom = (real)1 / (va + vb + vc); real vv = vb * denom, ww = vc * denom;
    u = 1 - vv - ww; v = vv; w = ww; return a + ab * vv + ac * ww;
}

inline bool sameSideFace(const Vector3& a, const Vector3& b, const Vector3& c,
                         const Vector3& ref, const Vector3& p) {
    Vector3 n = (b - a) % (c - a);
    return (n * (ref - a)) * (n * (p - a)) >= 0;      // p on same side of abc as ref
}
inline bool originInsideTetra(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    Vector3 O;
    return sameSideFace(a, b, c, d, O) && sameSideFace(a, c, d, b, O)
        && sameSideFace(a, d, b, c, O) && sameSideFace(b, c, d, a, O);
}

// Keep the simplex vertices whose barycentric weight is non-trivial.
inline void reduceKeep3(std::vector<Support>& Q, real w[4], real u, real v, real ww) {
    Support q0 = Q[0], q1 = Q[1], q2 = Q[2];
    std::vector<Support> nq; real nw[4] = {0, 0, 0, 0}; int n = 0;
    const real eps = (real)1e-9;
    if (u  > eps) { nq.push_back(q0); nw[n++] = u; }
    if (v  > eps) { nq.push_back(q1); nw[n++] = v; }
    if (ww > eps) { nq.push_back(q2); nw[n++] = ww; }
    if (n == 0) { nq.push_back(q0); nw[0] = 1; n = 1; }
    Q = nq; for (int i = 0; i < 4; i++) w[i] = i < n ? nw[i] : 0;
}

// Closest point on the current simplex (1..4 pts) to the origin; reduces the
// simplex to the supporting sub-feature and fills barycentric weights. Sets
// `inside` when a tetrahedron already encloses the origin (→ intersection).
inline Vector3 closestSimplex(std::vector<Support>& Q, real w[4], bool& inside) {
    inside = false; w[0] = w[1] = w[2] = w[3] = 0;
    if (Q.size() == 1) { w[0] = 1; return Q[0].v; }
    if (Q.size() == 2) {
        Vector3 a = Q[0].v, b = Q[1].v, ab = b - a; real den = ab * ab;
        real t = den > (real)1e-18 ? (-(a * ab)) / den : 0;
        if (t <= 0) { Q = {Q[0]}; w[0] = 1; return a; }
        if (t >= 1) { Q = {Q[1]}; w[0] = 1; return b; }
        w[0] = 1 - t; w[1] = t; return a + ab * t;
    }
    if (Q.size() == 3) {
        real u, v, ww; Vector3 c = closestPtTriangle(Vector3(), Q[0].v, Q[1].v, Q[2].v, u, v, ww);
        reduceKeep3(Q, w, u, v, ww); return c;
    }
    // 4 points
    if (originInsideTetra(Q[0].v, Q[1].v, Q[2].v, Q[3].v)) { inside = true; w[0] = 1; return Vector3(); }
    const int fc[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    real bestSq = REAL_MAX; Vector3 bestC; int bf = 0; real bu = 0, bv = 0, bw = 0;
    for (int f = 0; f < 4; f++) {
        real u, v, ww; Vector3 c = closestPtTriangle(Vector3(), Q[fc[f][0]].v, Q[fc[f][1]].v, Q[fc[f][2]].v, u, v, ww);
        real sq = c.squareMagnitude();
        if (sq < bestSq) { bestSq = sq; bestC = c; bf = f; bu = u; bv = v; bw = ww; }
    }
    std::vector<Support> face = {Q[fc[bf][0]], Q[fc[bf][1]], Q[fc[bf][2]]};
    Q = face; reduceKeep3(Q, w, bu, bv, bw); return bestC;
}

inline bool sameDir(const Vector3& a, const Vector3& b) { return a * b > 0; }

// Any non-zero vector perpendicular to v (cross with the least-aligned axis).
inline Vector3 anyPerp(const Vector3& v) {
    real ax = real_abs(v.x), ay = real_abs(v.y), az = real_abs(v.z);
    Vector3 t = (ax <= ay && ax <= az) ? Vector3(1, 0, 0) : (ay <= az ? Vector3(0, 1, 0) : Vector3(0, 0, 1));
    Vector3 p = v % t;
    if (p.squareMagnitude() < (real)1e-18) p = v % Vector3(0, 1, 0);
    return p;
}
// Direction along `edge`'s plane pointing toward the origin; a perpendicular
// when the origin lies on the edge's line (keeps GJK from collapsing to 0).
inline Vector3 towardOrigin(const Vector3& edge, const Vector3& ao) {
    Vector3 c = edge % ao;
    if (c.squareMagnitude() < (real)1e-18) return anyPerp(edge);
    return c % edge;
}

// --- GJK boolean: evolve a simplex toward the origin; on `true` the simplex is
// a tetrahedron enclosing the origin (fed straight to EPA). Muratori-style, with
// degeneracy guards so axis-aligned / collinear cases still reach a tetrahedron.
inline bool doSimplex(std::vector<Support>& s, Vector3& dir) {
    if (s.size() == 2) {
        Vector3 a = s[0].v, b = s[1].v, ab = b - a, ao = a * (real)-1;
        if (sameDir(ab, ao)) dir = towardOrigin(ab, ao); else { s = {s[0]}; dir = ao; }
        return false;
    }
    if (s.size() == 3) {
        Vector3 a = s[0].v, b = s[1].v, c = s[2].v;
        Vector3 ab = b - a, ac = c - a, ao = a * (real)-1, abc = ab % ac;
        if (abc.squareMagnitude() < (real)1e-18) { s = {s[0], s[1]}; return doSimplex(s, dir); } // collinear
        if (sameDir(abc % ac, ao)) {
            if (sameDir(ac, ao)) { s = {s[0], s[2]}; dir = towardOrigin(ac, ao); }
            else { s = {s[0], s[1]}; return doSimplex(s, dir); }
        } else if (sameDir(ab % abc, ao)) { s = {s[0], s[1]}; return doSimplex(s, dir); }
        else if (sameDir(abc, ao)) dir = abc;
        else { s = {s[0], s[2], s[1]}; dir = abc * (real)-1; }
        return false;
    }
    // tetrahedron
    Vector3 a = s[0].v, b = s[1].v, c = s[2].v, d = s[3].v;
    Vector3 ab = b - a, ac = c - a, ad = d - a, ao = a * (real)-1;
    Vector3 abc = ab % ac, acd = ac % ad, adb = ad % ab;
    if (sameDir(abc, ao)) { s = {s[0], s[1], s[2]}; return doSimplex(s, dir); }
    if (sameDir(acd, ao)) { s = {s[0], s[2], s[3]}; return doSimplex(s, dir); }
    if (sameDir(adb, ao)) { s = {s[0], s[3], s[1]}; return doSimplex(s, dir); }
    return true;                                                  // origin enclosed
}

inline bool gjkIntersect(const CollisionConvex& A, const CollisionConvex& B, std::vector<Support>& simplex) {
    Vector3 dir(1, 0, 0);
    Support s = mdSupport(A, B, dir);
    simplex.clear(); simplex.push_back(s);
    dir = s.v * (real)-1;
    for (int iter = 0; iter < 64; iter++) {
        if (dir.squareMagnitude() < (real)1e-14) return true;    // origin on the simplex (touching)
        s = mdSupport(A, B, dir);
        if (s.v * dir < 0) return false;                         // support didn't pass the origin
        simplex.insert(simplex.begin(), s);
        if (doSimplex(simplex, dir)) return true;
    }
    return false;
}

// --- GJK distance: minimum separation of two disjoint hulls (0 if intersecting),
// plus the witness (closest) points on each hull.
inline real gjkDistance(const CollisionConvex& A, const CollisionConvex& B,
                        Vector3* outA = nullptr, Vector3* outB = nullptr) {
    std::vector<Support> Q; Q.push_back(mdSupport(A, B, Vector3(1, 0, 0)));
    real w[4] = {1, 0, 0, 0};
    Vector3 v = Q[0].v;
    for (int iter = 0; iter < 64; iter++) {
        Vector3 dir = v * (real)-1;
        if (v.squareMagnitude() < (real)1e-18) { v.clear(); break; }   // origin reached → 0
        Support s = mdSupport(A, B, dir);
        real vv = v * v, vw = v * s.v;
        if (vv - vw <= (real)1e-12 * (vv + 1)) break;                  // support gains nothing → converged
        Q.insert(Q.begin(), s);
        bool inside = false; v = closestSimplex(Q, w, inside);
        if (inside) { v.clear(); break; }
    }
    if (outA || outB) {
        Vector3 cA, cB;
        for (size_t i = 0; i < Q.size(); i++) { cA += Q[i].a * w[i]; cB += Q[i].b * w[i]; }
        if (outA) *outA = cA;
        if (outB) *outB = cB;
    }
    return v.magnitude();
}

// --- EPA: expand the GJK tetrahedron to the closest face of the Minkowski
// boundary → penetration depth, outward normal, and a witness contact point.
inline bool epa(const CollisionConvex& A, const CollisionConvex& B, std::vector<Support>& verts,
                Vector3& outNormal, real& outDepth, Vector3& outPoint) {
    struct Face { int a, b, c; Vector3 n; real d; };
    std::vector<Face> faces;
    // Orient every face outward from a fixed interior point (the initial-tetra
    // centroid, strictly inside and staying inside as the polytope only grows).
    // Orienting by the origin instead breaks when the origin lies on a face —
    // exactly the degenerate case of two aligned boxes.
    Vector3 interior = (verts[0].v + verts[1].v + verts[2].v + verts[3].v) * (real)0.25;
    auto addFace = [&](int i, int j, int k) {
        Vector3 n = (verts[j].v - verts[i].v) % (verts[k].v - verts[i].v);
        real len = n.magnitude(); if (len < (real)1e-12) return;
        n = n * ((real)1 / len);
        if (n * (verts[i].v - interior) < 0) { n = n * (real)-1; int t = j; j = k; k = t; }
        faces.push_back({i, j, k, n, n * verts[i].v});      // d = signed origin-distance
    };
    addFace(0, 1, 2); addFace(0, 1, 3); addFace(0, 2, 3); addFace(1, 2, 3);
    if (faces.size() < 4) return false;

    Face best{};
    for (int iter = 0; iter < 64; iter++) {
        int bi = 0; for (int i = 1; i < (int)faces.size(); i++) if (faces[i].d < faces[bi].d) bi = i;
        best = faces[bi];
        Support s = mdSupport(A, B, best.n);
        if (s.v * best.n - best.d < (real)1e-4) break;             // can't push the face out → done
        int ni = (int)verts.size(); verts.push_back(s);
        std::vector<std::pair<int, int>> edges;
        auto addEdge = [&](int a, int b) {
            for (size_t e = 0; e < edges.size(); e++)
                if (edges[e].first == b && edges[e].second == a) { edges.erase(edges.begin() + e); return; }
            edges.push_back({a, b});
        };
        for (int i = (int)faces.size() - 1; i >= 0; i--)
            if (faces[i].n * (s.v - verts[faces[i].a].v) > (real)1e-9) {   // face visible from new point
                addEdge(faces[i].a, faces[i].b); addEdge(faces[i].b, faces[i].c); addEdge(faces[i].c, faces[i].a);
                faces.erase(faces.begin() + i);
            }
        for (auto& e : edges) addFace(e.first, e.second, ni);
        if (faces.empty()) return false;
    }
    outNormal = best.n; outDepth = best.d;
    Vector3 q = best.n * best.d;                                    // origin projected onto the face
    real u, v, w; closestPtTriangle(q, verts[best.a].v, verts[best.b].v, verts[best.c].v, u, v, w);
    Vector3 cA = verts[best.a].a * u + verts[best.b].a * v + verts[best.c].a * w;
    Vector3 cB = verts[best.a].b * u + verts[best.b].b * v + verts[best.c].b * w;
    outPoint = (cA + cB) * (real)0.5;
    return true;
}

inline void perpBasis(const Vector3& ax, Vector3& u, Vector3& v) {
    Vector3 t = real_abs(ax.x) < (real)0.9 ? Vector3(1, 0, 0) : Vector3(0, 1, 0);
    u = ax % t; u.normalise(); v = ax % u; v.normalise();
}

} // namespace cvx

// ------------------------------------------------------------- static tri-mesh
// A triangle soup used as an immovable collider (terrain, level geometry). Stores
// world-space triangles with precomputed normals + AABBs; build() adds a uniform
// grid so queries touch only nearby triangles instead of the whole soup.
class CollisionTriMesh {
public:
    struct Tri { Vector3 a, b, c, normal, lo, hi; };
    std::vector<Tri> tris;
    Vector3 lo, hi;                                                // overall AABB
    bool haveBounds = false;

    void addTriangle(const Vector3& a, const Vector3& b, const Vector3& c) {
        Tri t; t.a = a; t.b = b; t.c = c;
        t.normal = (b - a) % (c - a); t.normal.normalise();
        for (int i = 0; i < 3; i++) {
            real mn = a[i], mx = a[i];
            mn = b[i] < mn ? b[i] : mn; mx = b[i] > mx ? b[i] : mx;
            mn = c[i] < mn ? c[i] : mn; mx = c[i] > mx ? c[i] : mx;
            t.lo[i] = mn; t.hi[i] = mx;
        }
        if (!haveBounds) { lo = t.lo; hi = t.hi; haveBounds = true; }
        else for (int i = 0; i < 3; i++) { if (t.lo[i] < lo[i]) lo[i] = t.lo[i]; if (t.hi[i] > hi[i]) hi[i] = t.hi[i]; }
        tris.push_back(t);
    }

    // Build the optional uniform-grid midphase (~targetDim cells across the longest axis).
    void build(int targetDim = 16) {
        cells.clear();
        Vector3 ext = hi - lo;
        real m = ext.x; if (ext.y > m) m = ext.y; if (ext.z > m) m = ext.z;
        cell = m > 0 ? m / targetDim : 1; if (cell <= 0) cell = 1;
        for (int i = 0; i < (int)tris.size(); i++) {
            int x0 = cc(tris[i].lo.x - lo.x), x1 = cc(tris[i].hi.x - lo.x);
            int y0 = cc(tris[i].lo.y - lo.y), y1 = cc(tris[i].hi.y - lo.y);
            int z0 = cc(tris[i].lo.z - lo.z), z1 = cc(tris[i].hi.z - lo.z);
            for (int x = x0; x <= x1; x++) for (int y = y0; y <= y1; y++) for (int z = z0; z <= z1; z++)
                cells[key(x, y, z)].push_back(i);
        }
    }

    // Candidate triangle indices whose cells overlap [qlo,qhi]. Falls back to the
    // full list when no grid was built. Never misses a real overlap (extra
    // candidates are cheap — the detector AABB-rejects them).
    void queryCandidates(const Vector3& qlo, const Vector3& qhi, std::vector<int>& out) const {
        out.clear();
        if (cells.empty()) { for (int i = 0; i < (int)tris.size(); i++) out.push_back(i); return; }
        int x0 = cc(qlo.x - lo.x), x1 = cc(qhi.x - lo.x);
        int y0 = cc(qlo.y - lo.y), y1 = cc(qhi.y - lo.y);
        int z0 = cc(qlo.z - lo.z), z1 = cc(qhi.z - lo.z);
        std::vector<char> seen(tris.size(), 0);
        for (int x = x0; x <= x1; x++) for (int y = y0; y <= y1; y++) for (int z = z0; z <= z1; z++) {
            auto it = cells.find(key(x, y, z)); if (it == cells.end()) continue;
            for (int idx : it->second) if (!seen[idx]) { seen[idx] = 1; out.push_back(idx); }
        }
    }
private:
    std::unordered_map<long long, std::vector<int>> cells;
    real cell = 1;
    int cc(real rel) const { return (int)std::floor(rel / cell); }
    static long long key(int x, int y, int z) {
        auto clampf = [](int v) { if (v < -32000) v = -32000; if (v > 32000) v = 32000; return (long long)(v + 32768); };
        return (clampf(x) << 32) | (clampf(y) << 16) | clampf(z);
    }
};

// ------------------------------------------------------------- detectors
struct ConvexCollision {
    // ---- convex ↔ convex (GJK + EPA) ----
    static bool intersect(const CollisionConvex& a, const CollisionConvex& b) {
        std::vector<cvx::Support> s; return cvx::gjkIntersect(a, b, s);
    }
    static real distance(const CollisionConvex& a, const CollisionConvex& b,
                         Vector3* closestA = nullptr, Vector3* closestB = nullptr) {
        return cvx::gjkDistance(a, b, closestA, closestB);
    }
    static unsigned convexAndConvex(const CollisionConvex& a, const CollisionConvex& b, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        std::vector<cvx::Support> simplex;
        if (!cvx::gjkIntersect(a, b, simplex) || simplex.size() < 4) return 0;
        Vector3 n, point; real depth;
        if (!cvx::epa(a, b, simplex, n, depth, point) || depth <= 0) return 0;
        Contact* c = data->contacts;
        c->contactNormal = n * (real)-1;              // EPA normal is A→out; flip to point body[1]→body[0]
        c->penetration = depth;
        c->contactPoint = point;
        c->setBodyData(a.body, b.body, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }

    // ---- shared approximations for round primitives ----
    // Treat a rounded segment (capsule of the given radius) against a sphere/box.
    static unsigned segSphere(const Vector3& a, const Vector3& b, real R, RigidBody* body,
                              const CollisionSphere& sph, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 sc = sph.getAxis(3);
        Vector3 p = cvx::closestOnSeg(a, b, sc);
        Vector3 mid = p - sc; real dist = mid.magnitude();
        if (dist <= 0 || dist >= R + sph.radius) return 0;
        Contact* c = data->contacts;
        c->contactNormal = mid * ((real)1 / dist);    // sphere(body[1]) → primitive(body[0])
        c->penetration = R + sph.radius - dist;
        c->contactPoint = sc + mid * (sph.radius / (R + sph.radius));
        c->setBodyData(body, sph.body, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    static unsigned segBox(const Vector3& a, const Vector3& b, real R, RigidBody* body,
                           const CollisionBox& box, CollisionData* data) {
        unsigned used = 0;
        for (int pass = 0; pass < 2 && data->hasMoreContacts(); pass++) {
            Vector3 probe = pass == 0 ? a : b;
            Vector3 rel = box.getTransform().transformInverse(probe), clamped = rel;
            for (int i = 0; i < 3; i++) { if (clamped[i] > box.halfSize[i]) clamped[i] = box.halfSize[i];
                if (clamped[i] < -box.halfSize[i]) clamped[i] = -box.halfSize[i]; }
            Vector3 cw = box.getTransform().transform(clamped);
            Vector3 seg = cvx::closestOnSeg(a, b, cw);
            Vector3 d = seg - cw; real dist = d.magnitude();
            if (dist <= 0 || dist >= R) continue;
            Contact* c = data->contacts;
            c->contactNormal = d * ((real)1 / dist);
            c->penetration = R - dist;
            c->contactPoint = cw;
            c->setBodyData(body, box.body, data->friction, data->restitution);
            data->addContacts(1); used++;
        }
        return used;
    }

    // ---- cylinder ↔ plane / sphere / box ----
    static unsigned cylinderAndHalfSpace(const CollisionCylinder& cyl, const CollisionPlane& plane, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 n = plane.direction, c = cyl.getAxis(3), ax = cyl.getAxis(1), u, v;
        cvx::perpBasis(ax, u, v);
        unsigned used = 0;
        auto tryPt = [&](const Vector3& p) {
            if (!data->hasMoreContacts()) return;
            real d = n * p - plane.offset; if (d >= 0) return;
            Contact* ct = data->contacts;
            ct->contactNormal = n; ct->penetration = -d; ct->contactPoint = p;
            ct->setBodyData(cyl.body, nullptr, data->friction, data->restitution);
            data->addContacts(1); used++;
        };
        for (int cap = 0; cap < 2; cap++) {
            Vector3 cc = c + ax * (cap ? cyl.halfHeight : -cyl.halfHeight);
            for (int k = 0; k < 8; k++) { real ang = 2 * real_pi * k / 8;
                tryPt(cc + (u * real_cos(ang) + v * real_sin(ang)) * cyl.radius); }
        }
        tryPt(cyl.support(n * (real)-1));                 // deepest point (covers side rests)
        return used;
    }
    static unsigned cylinderAndSphere(const CollisionCylinder& cyl, const CollisionSphere& sph, CollisionData* data) {
        Vector3 a, b; cyl.endpoints(a, b);
        return segSphere(a, b, cyl.radius, cyl.body, sph, data);
    }
    static unsigned cylinderAndBox(const CollisionCylinder& cyl, const CollisionBox& box, CollisionData* data) {
        Vector3 a, b; cyl.endpoints(a, b);
        return segBox(a, b, cyl.radius, cyl.body, box, data);
    }

    // ---- cone ↔ plane / sphere / box ----
    static unsigned coneAndHalfSpace(const CollisionCone& cone, const CollisionPlane& plane, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 n = plane.direction, ax = cone.getAxis(1), u, v;
        cvx::perpBasis(ax, u, v);
        Vector3 base = cone.baseCentre();
        unsigned used = 0;
        auto tryPt = [&](const Vector3& p) {
            if (!data->hasMoreContacts()) return;
            real d = n * p - plane.offset; if (d >= 0) return;
            Contact* ct = data->contacts;
            ct->contactNormal = n; ct->penetration = -d; ct->contactPoint = p;
            ct->setBodyData(cone.body, nullptr, data->friction, data->restitution);
            data->addContacts(1); used++;
        };
        for (int k = 0; k < 8; k++) { real ang = 2 * real_pi * k / 8;
            tryPt(base + (u * real_cos(ang) + v * real_sin(ang)) * cone.radius); }
        tryPt(cone.apex());
        tryPt(cone.support(n * (real)-1));
        return used;
    }
    // cone approximated by a rounded axis segment (base radius) for sphere/box.
    static unsigned coneAndSphere(const CollisionCone& cone, const CollisionSphere& sph, CollisionData* data) {
        return segSphere(cone.apex(), cone.baseCentre(), cone.radius, cone.body, sph, data);
    }
    static unsigned coneAndBox(const CollisionCone& cone, const CollisionBox& box, CollisionData* data) {
        return segBox(cone.apex(), cone.baseCentre(), cone.radius, cone.body, box, data);
    }

    // ---- static tri-mesh ↔ sphere / box ----
    static unsigned sphereAndTriMesh(const CollisionSphere& sph, const CollisionTriMesh& mesh, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 sc = sph.getAxis(3); real r = sph.radius;
        Vector3 qlo = sc - Vector3(r, r, r), qhi = sc + Vector3(r, r, r);
        std::vector<int> cand; mesh.queryCandidates(qlo, qhi, cand);
        real bestDist = REAL_MAX; int bestTri = -1; Vector3 bestP;
        for (int idx : cand) {
            const CollisionTriMesh::Tri& t = mesh.tris[idx];
            if (qhi.x < t.lo.x || qlo.x > t.hi.x || qhi.y < t.lo.y ||
                qlo.y > t.hi.y || qhi.z < t.lo.z || qlo.z > t.hi.z) continue;
            real u, v, w; Vector3 p = cvx::closestPtTriangle(sc, t.a, t.b, t.c, u, v, w);
            real d = (sc - p).magnitude();
            if (d < r && d < bestDist) { bestDist = d; bestTri = idx; bestP = p; }
        }
        if (bestTri < 0) return 0;                          // closest triangle → single contact
        const CollisionTriMesh::Tri& t = mesh.tris[bestTri];
        Vector3 diff = sc - bestP; real dist = diff.magnitude();
        Contact* c = data->contacts;
        c->contactNormal = dist > (real)1e-9 ? diff * ((real)1 / dist) : t.normal;  // surface → sphere
        c->penetration = r - dist;
        c->contactPoint = bestP;
        c->setBodyData(sph.body, nullptr, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    static unsigned boxAndTriMesh(const CollisionBox& box, const CollisionTriMesh& mesh, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 c = box.getAxis(3);
        Vector3 ext(transformToAxis(box, Vector3(1, 0, 0)),
                    transformToAxis(box, Vector3(0, 1, 0)),
                    transformToAxis(box, Vector3(0, 0, 1)));
        Vector3 qlo = c - ext, qhi = c + ext;
        std::vector<int> cand; mesh.queryCandidates(qlo, qhi, cand);
        static const real mult[8][3] = {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
        unsigned used = 0;
        for (int idx : cand) {
            if (!data->hasMoreContacts()) break;
            const CollisionTriMesh::Tri& t = mesh.tris[idx];
            if (qhi.x < t.lo.x || qlo.x > t.hi.x || qhi.y < t.lo.y ||
                qlo.y > t.hi.y || qhi.z < t.lo.z || qlo.z > t.hi.z) continue;
            Vector3 N = t.normal; if (N * (c - t.a) < 0) N = N * (real)-1;   // orient toward the box
            real planeOff = N * t.a;
            for (int i = 0; i < 8 && data->hasMoreContacts(); i++) {
                Vector3 v(mult[i][0] * box.halfSize.x, mult[i][1] * box.halfSize.y, mult[i][2] * box.halfSize.z);
                Vector3 corner = box.getTransform().transform(v);
                real sd = N * corner - planeOff; if (sd >= 0) continue;      // above the triangle plane
                Vector3 proj = corner - N * sd;                             // corner lifted onto the plane
                real u, vv, w; Vector3 cp = cvx::closestPtTriangle(proj, t.a, t.b, t.c, u, vv, w);
                if ((cp - proj).squareMagnitude() > (real)1e-8) continue;    // projects outside the triangle
                Contact* ct = data->contacts;
                ct->contactNormal = N; ct->penetration = -sd; ct->contactPoint = proj;
                ct->setBodyData(box.body, nullptr, data->friction, data->restitution);
                data->addContacts(1); used++;
            }
        }
        return used;
    }
};

} // namespace phys
