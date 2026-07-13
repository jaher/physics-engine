// Scene queries: raycasts against sphere / box / plane / capsule primitives, and
// sphere-overlap tests. The query API every surveyed engine exposes.
#pragma once
#include "shapes.h"

namespace phys {

struct Ray { Vector3 origin, dir; };                      // dir need not be unit (normalised on use)
struct RayHit { real t = REAL_MAX; Vector3 point, normal; RigidBody* body = nullptr; bool hit = false; };

inline bool raySphere(const Ray& r, const Vector3& c, real rad, RayHit& out) {
    Vector3 d = r.dir.unit(), m = r.origin - c;
    real b = m * d, cc = m * m - rad * rad;
    if (cc > 0 && b > 0) return false;
    real disc = b * b - cc; if (disc < 0) return false;
    real t = -b - real_sqrt(disc); if (t < 0) t = 0;
    out.hit = true; out.t = t; out.point = r.origin + d * t;
    out.normal = (out.point - c) * (((real)1) / rad);
    return true;
}
inline bool rayPlane(const Ray& r, const Vector3& n, real offset, RayHit& out) {
    Vector3 d = r.dir.unit();
    real denom = n * d; if (real_abs(denom) < 1e-12) return false;
    real t = (offset - n * r.origin) / denom; if (t < 0) return false;
    out.hit = true; out.t = t; out.point = r.origin + d * t;
    out.normal = denom < 0 ? n : n * -1;
    return true;
}
inline bool rayBox(const Ray& r, const CollisionBox& box, RayHit& out) {
    Vector3 d = r.dir.unit();
    Vector3 o = box.getTransform().transformInverse(r.origin);
    Vector3 dl = box.getTransform().transformInverseDirection(d);
    real tmin = 0, tmax = REAL_MAX; int axisMin = -1; real signMin = 1;
    for (int i = 0; i < 3; i++) {
        if (real_abs(dl[i]) < 1e-12) { if (o[i] < -box.halfSize[i] || o[i] > box.halfSize[i]) return false; continue; }
        real inv = ((real)1) / dl[i];
        real t1 = (-box.halfSize[i] - o[i]) * inv, t2 = (box.halfSize[i] - o[i]) * inv;
        real sgn = -1; if (t1 > t2) { real tmp = t1; t1 = t2; t2 = tmp; sgn = 1; }
        if (t1 > tmin) { tmin = t1; axisMin = i; signMin = sgn; }
        if (t2 < tmax) tmax = t2;
        if (tmin > tmax) return false;
    }
    if (axisMin < 0) return false;
    out.hit = true; out.t = tmin; out.point = r.origin + d * tmin;
    Vector3 nl; nl[axisMin] = signMin;
    out.normal = box.getTransform().transformDirection(nl);
    return true;
}
inline bool rayCapsule(const Ray& r, const CollisionCapsule& cap, RayHit& out) {
    // conservative: march the ray against distance-to-segment (sphere-traced)
    Vector3 a, b; cap.endpoints(a, b);
    Vector3 d = r.dir.unit();
    real t = 0;
    for (int i = 0; i < 128; i++) {
        Vector3 p = r.origin + d * t;
        Vector3 q = closestOnSegment(a, b, p);
        real dist = (p - q).magnitude() - cap.radius;
        if (dist < 1e-5) { out.hit = true; out.t = t; out.point = p;
            out.normal = (p - q).unit(); return true; }
        t += dist; if (t > 1e6) break;
    }
    return false;
}

// Aggregate query over a set of primitives (each may carry a body pointer).
struct RaycastWorld {
    std::vector<CollisionSphere*> spheres;
    std::vector<CollisionBox*> boxes;
    std::vector<CollisionCapsule*> capsules;
    std::vector<CollisionPlane*> planes;
    RayHit raycast(const Ray& r) const {
        RayHit best;
        RayHit h;
        for (auto* s : spheres) if (raySphere(r, s->getAxis(3), s->radius, h) && h.t < best.t) { best = h; best.body = s->body; }
        for (auto* b : boxes) if (rayBox(r, *b, h) && h.t < best.t) { best = h; best.body = b->body; }
        for (auto* c : capsules) if (rayCapsule(r, *c, h) && h.t < best.t) { best = h; best.body = c->body; }
        for (auto* p : planes) if (rayPlane(r, p->direction, p->offset, h) && h.t < best.t) { best = h; best.body = nullptr; }
        return best;
    }
};

} // namespace phys
