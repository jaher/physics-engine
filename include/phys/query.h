// Scene-query header: shape *sweep* queries (swept sphere / swept AABB → first
// time-of-impact), *overlap* queries (all bodies whose bounds intersect a query
// sphere / AABB) and a *speculative-contact* helper (predict a contact for a
// fast body before it penetrates, so the solver can stop it this step). These are
// the scene-query APIs every mainstream engine exposes (PhysX sweep/overlap,
// Bullet contactTest/convexSweepTest, Unity Physics.SphereCast/OverlapBox).
#pragma once
#include "shapes.h"          // CollisionSphere/Box/Plane/Capsule + collide_fine + contacts
#include <algorithm>
#include <vector>

namespace phys {

// ---------------------------------------------------------------- axis-aligned box
struct Aabb {
    Vector3 min, max;
    Aabb() : min(REAL_MAX, REAL_MAX, REAL_MAX), max(-REAL_MAX, -REAL_MAX, -REAL_MAX) {}
    Aabb(const Vector3& mn, const Vector3& mx) : min(mn), max(mx) {}
    static Aabb fromCentreHalf(const Vector3& c, const Vector3& h) { return Aabb(c - h, c + h); }

    Vector3 centre()  const { return (min + max) * (real)0.5; }
    Vector3 extents() const { return (max - min) * (real)0.5; }

    // inclusive overlap (touching counts) — same predicate SAP & the tree share
    bool overlaps(const Aabb& o) const {
        return min.x <= o.max.x && max.x >= o.min.x &&
               min.y <= o.max.y && max.y >= o.min.y &&
               min.z <= o.max.z && max.z >= o.min.z;
    }
    bool contains(const Aabb& o) const {
        return min.x <= o.min.x && max.x >= o.max.x &&
               min.y <= o.min.y && max.y >= o.max.y &&
               min.z <= o.min.z && max.z >= o.max.z;
    }
    Aabb merged(const Aabb& o) const {
        return Aabb(Vector3(std::min(min.x, o.min.x), std::min(min.y, o.min.y), std::min(min.z, o.min.z)),
                    Vector3(std::max(max.x, o.max.x), std::max(max.y, o.max.y), std::max(max.z, o.max.z)));
    }
    Aabb expanded(real margin) const {
        Vector3 m(margin, margin, margin);
        return Aabb(min - m, max + m);
    }
    // surface area (SAH cost metric for the dynamic tree) and volume
    real area()   const { Vector3 d = max - min; return 2 * (d.x * d.y + d.y * d.z + d.z * d.x); }
    real volume() const { Vector3 d = max - min; return d.x * d.y * d.z; }

    // exterior distance from a point to the box (0 if inside)
    real distanceTo(const Vector3& p) const {
        real dx = std::max((real)0, std::max(min.x - p.x, p.x - max.x));
        real dy = std::max((real)0, std::max(min.y - p.y, p.y - max.y));
        real dz = std::max((real)0, std::max(min.z - p.z, p.z - max.z));
        return real_sqrt(dx * dx + dy * dy + dz * dz);
    }
    bool overlapsSphere(const Vector3& c, real r) const { return distanceTo(c) <= r; }
};

// world-space AABB of the collision primitives ------------------------------------
inline Aabb aabbOfSphere(const CollisionSphere& s) {
    Vector3 c = s.getAxis(3), r(s.radius, s.radius, s.radius);
    return Aabb(c - r, c + r);
}
inline Aabb aabbOfBox(const CollisionBox& b) {
    Vector3 c = b.getAxis(3);
    // extent along each world axis = Σ |R_row · e_k| * halfSize   (|R|·h)
    Vector3 ext(
        real_abs(b.getAxis(0).x) * b.halfSize.x + real_abs(b.getAxis(1).x) * b.halfSize.y + real_abs(b.getAxis(2).x) * b.halfSize.z,
        real_abs(b.getAxis(0).y) * b.halfSize.x + real_abs(b.getAxis(1).y) * b.halfSize.y + real_abs(b.getAxis(2).y) * b.halfSize.z,
        real_abs(b.getAxis(0).z) * b.halfSize.x + real_abs(b.getAxis(1).z) * b.halfSize.y + real_abs(b.getAxis(2).z) * b.halfSize.z);
    return Aabb(c - ext, c + ext);
}

// ------------------------------------------------------------------ sweep results
struct SweepHit {
    bool hit = false;
    real toi = REAL_MAX;         // fraction of the sweep in [0,1] at first contact
    Vector3 point, normal;       // world contact point + surface normal (toward the mover)
    RigidBody* body = nullptr;   // the shape that was hit (nullptr for planes)
};

// swept sphere vs an (infinite) half-space plane — analytic, exact -----------------
inline bool sweptSphereVsPlane(const Vector3& p0, const Vector3& p1, real r,
                               const Vector3& n, real offset, real& toi, Vector3& hitN) {
    real d0 = n * p0 - offset;             // signed distance of the centre at t=0
    real d1 = n * p1 - offset;             //                          ... at t=1
    if (d0 <= r) { toi = 0; hitN = n; return true; }   // already within a radius
    if (d1 > r)  return false;             // never reaches within a radius
    toi = (d0 - r) / (d0 - d1);            // linear crossing of the d==r shell
    hitN = n;
    return true;
}
// swept sphere vs a static sphere — analytic (ray vs Minkowski-summed sphere) -------
inline bool sweptSphereVsSphere(const Vector3& p0, const Vector3& p1, real r,
                                const Vector3& c, real cr, real& toi, Vector3& hitN) {
    Vector3 d = p1 - p0, m = p0 - c; real R = r + cr;
    real a = d * d;
    real cc = m * m - R * R;
    if (cc <= 0) { toi = 0; hitN = (p0 - c).unit(); return true; }   // start overlapping
    if (a < (real)1e-18) return false;                               // not moving
    real b = m * d;
    if (b > 0) return false;                                         // moving away
    real disc = b * b - a * cc;
    if (disc < 0) return false;
    real tt = (-b - real_sqrt(disc)) / a;
    if (tt < 0 || tt > 1) return false;
    toi = tt; hitN = ((p0 + d * tt) - c).unit();
    return true;
}
// swept sphere vs an oriented box — conservative advancement (sphere-cast). Handles
// rounded corners exactly and never tunnels: marches by the exact exterior gap.
inline bool sweptSphereVsBox(const Vector3& p0, const Vector3& p1, real r,
                             const CollisionBox& box, real& toi, Vector3& hitN) {
    Vector3 seg = p1 - p0; real L = seg.magnitude();
    const Matrix4& xf = box.getTransform();
    auto gap = [&](const Vector3& worldC) -> real {   // exterior distance centre→box − r
        Vector3 rel = xf.transformInverse(worldC);
        Vector3 d(std::max((real)0, real_abs(rel.x) - box.halfSize.x),
                  std::max((real)0, real_abs(rel.y) - box.halfSize.y),
                  std::max((real)0, real_abs(rel.z) - box.halfSize.z));
        return d.magnitude() - r;
    };
    auto normalAt = [&](const Vector3& worldC) -> Vector3 {
        Vector3 rel = xf.transformInverse(worldC), clamped = rel;
        for (int i = 0; i < 3; i++) { if (clamped[i] >  box.halfSize[i]) clamped[i] =  box.halfSize[i];
                                      if (clamped[i] < -box.halfSize[i]) clamped[i] = -box.halfSize[i]; }
        Vector3 nl = (rel - clamped);
        if (nl.squareMagnitude() < (real)1e-18) {           // centre inside the box: push out on nearest face
            real bx = box.halfSize.x - real_abs(rel.x), by = box.halfSize.y - real_abs(rel.y), bz = box.halfSize.z - real_abs(rel.z);
            if (bx <= by && bx <= bz) nl = Vector3(rel.x < 0 ? -1 : 1, 0, 0);
            else if (by <= bz)        nl = Vector3(0, rel.y < 0 ? -1 : 1, 0);
            else                      nl = Vector3(0, 0, rel.z < 0 ? -1 : 1);
        }
        return xf.transformDirection(nl).unit();
    };
    real g0 = gap(p0);
    if (g0 <= (real)1e-9) { toi = 0; hitN = normalAt(p0); return true; }   // already touching
    if (L < (real)1e-18) return false;                                     // not moving, and not touching
    real t = 0;
    for (int it = 0; it < 64; it++) {
        Vector3 c = p0 + seg * t;
        real g = gap(c);
        if (g <= (real)1e-7) { toi = t; hitN = normalAt(c); return true; }
        t += g / L;                        // safe step: cannot close more than g of clearance
        if (t > 1) return false;
    }
    return false;
}

// swept AABB vs a static AABB — slab method, first TOI + face normal ---------------
inline bool sweptAabbVsAabb(const Aabb& a, const Vector3& disp, const Aabb& b,
                            real& toi, Vector3& hitN) {
    real tEntry[3], tExit[3];
    for (int i = 0; i < 3; i++) {
        real amin = a.min[i], amax = a.max[i], bmin = b.min[i], bmax = b.max[i], v = disp[i];
        if (real_abs(v) < (real)1e-18) {
            if (amax < bmin || amin > bmax) return false;   // separated on this axis and not closing
            tEntry[i] = -REAL_MAX; tExit[i] = REAL_MAX;
        } else if (v > 0) {
            tEntry[i] = (bmin - amax) / v;
            tExit[i]  = (bmax - amin) / v;
        } else {
            tEntry[i] = (bmax - amin) / v;
            tExit[i]  = (bmin - amax) / v;
        }
    }
    real entry = std::max(tEntry[0], std::max(tEntry[1], tEntry[2]));
    real exitT = std::min(tExit[0],  std::min(tExit[1],  tExit[2]));
    if (entry > exitT) return false;
    if (entry < 0) {
        if (exitT <= 0) return false;   // already fully past
        toi = 0;                        // already overlapping at the start of the sweep
    } else {
        if (entry > 1) return false;    // impact is beyond this step
        toi = entry;
    }
    int axis = (tEntry[0] >= tEntry[1] && tEntry[0] >= tEntry[2]) ? 0 : (tEntry[1] >= tEntry[2] ? 1 : 2);
    hitN = Vector3();
    hitN[axis] = disp[axis] > 0 ? (real)-1 : (real)1;
    return true;
}

// -------------------------------------------------------------- speculative contact
// Predict a contact for a fast body *before* it penetrates. `separation` is the
// current gap (>0) or penetration (<0); `willHit` is true if the body's motion this
// step closes the gap. The resulting Contact lets the existing solver clamp the
// approach velocity (speculative contacts, à la Box2D / Bullet's CCD margins).
struct SpeculativeContact {
    bool willHit = false;
    real toi = 1;                 // fraction of dt at predicted impact
    real separation = 0;          // >0 gap, <=0 penetrating
    Vector3 point, normal;        // predicted surface point + normal toward the body
    RigidBody* body = nullptr;
    // Fold into a resolvable Contact: real overlap → true penetration; a mere gap →
    // zero penetration (the velocity solver still clamps the closing speed).
    Contact toContact(real restitution = 0, real friction = 0) const {
        Contact c;
        c.contactNormal = normal;
        c.contactPoint  = point;
        c.penetration   = separation < 0 ? -separation : (real)0;
        c.setBodyData(body, nullptr, friction, restitution);
        return c;
    }
};

inline SpeculativeContact speculativeSpherePlane(RigidBody* body, real radius,
                                                 const Vector3& n, real offset, real dt) {
    SpeculativeContact sc; sc.body = body; sc.normal = n;
    Vector3 p = body->getPosition();
    real dist = n * p - offset;                 // centre → plane (signed)
    sc.separation = dist - radius;              // gap to the sphere surface
    sc.point = p - n * radius;                  // near point of the sphere
    real vn = body->getVelocity() * n;          // normal velocity (<0 = approaching)
    if (sc.separation <= 0) { sc.willHit = true; sc.toi = 0; return sc; }
    real approach = -vn * dt;                   // distance closed toward the plane over dt
    if (approach > sc.separation) { sc.willHit = true; sc.toi = sc.separation / approach; }
    return sc;
}
inline SpeculativeContact speculativeSphereSphere(RigidBody* body, real radius,
                                                  const Vector3& staticC, real staticR, real dt) {
    SpeculativeContact sc; sc.body = body;
    Vector3 p = body->getPosition();
    Vector3 d = p - staticC; real dist = d.magnitude();
    Vector3 n = dist > 0 ? d * (((real)1) / dist) : Vector3(0, 1, 0);
    sc.normal = n;
    sc.separation = dist - (radius + staticR);
    sc.point = p - n * radius;
    real vn = body->getVelocity() * n;          // (<0 = approaching the static sphere)
    if (sc.separation <= 0) { sc.willHit = true; sc.toi = 0; return sc; }
    real approach = -vn * dt;
    if (approach > sc.separation) { sc.willHit = true; sc.toi = sc.separation / approach; }
    return sc;
}

// -------------------------------------------------------------------- QueryWorld
// A bag of static primitives supporting sweep and overlap queries (mirrors the
// existing RaycastWorld layout in raycast.h).
struct QueryWorld {
    std::vector<CollisionSphere*> spheres;
    std::vector<CollisionBox*>    boxes;
    std::vector<CollisionPlane*>  planes;    // planes participate in sphere-sweeps only

    // Sweep a sphere of radius `r` from `from` to `to`; earliest hit wins.
    SweepHit sweepSphere(const Vector3& from, const Vector3& to, real r) const {
        SweepHit best; real t; Vector3 nrm;
        for (auto* s : spheres)
            if (sweptSphereVsSphere(from, to, r, s->getAxis(3), s->radius, t, nrm) && t < best.toi) {
                best.hit = true; best.toi = t; best.normal = nrm; best.body = s->body; best.point = from + (to - from) * t + nrm * (-r);
            }
        for (auto* b : boxes)
            if (sweptSphereVsBox(from, to, r, *b, t, nrm) && t < best.toi) {
                best.hit = true; best.toi = t; best.normal = nrm; best.body = b->body; best.point = from + (to - from) * t + nrm * (-r);
            }
        for (auto* p : planes)
            if (sweptSphereVsPlane(from, to, r, p->direction, p->offset, t, nrm) && t < best.toi) {
                best.hit = true; best.toi = t; best.normal = nrm; best.body = nullptr; best.point = from + (to - from) * t + nrm * (-r);
            }
        return best;
    }

    // Sweep an AABB by displacement `disp` through the static shapes' bounds.
    SweepHit sweepAabb(const Aabb& box, const Vector3& disp) const {
        SweepHit best; real t; Vector3 nrm;
        for (auto* s : spheres)
            if (sweptAabbVsAabb(box, disp, aabbOfSphere(*s), t, nrm) && t < best.toi) {
                best.hit = true; best.toi = t; best.normal = nrm; best.body = s->body; best.point = box.centre() + disp * t;
            }
        for (auto* b : boxes)
            if (sweptAabbVsAabb(box, disp, aabbOfBox(*b), t, nrm) && t < best.toi) {
                best.hit = true; best.toi = t; best.normal = nrm; best.body = b->body; best.point = box.centre() + disp * t;
            }
        return best;
    }

    // Overlap: every body whose bounding shape intersects the query sphere / AABB.
    std::vector<RigidBody*> overlapSphere(const Vector3& c, real r) const {
        std::vector<RigidBody*> out;
        for (auto* s : spheres) {
            Vector3 m = s->getAxis(3) - c; real rr = s->radius + r;
            if (m.squareMagnitude() <= rr * rr) out.push_back(s->body);
        }
        for (auto* b : boxes)
            if (aabbOfBox(*b).overlapsSphere(c, r)) out.push_back(b->body);
        return out;
    }
    std::vector<RigidBody*> overlapAabb(const Aabb& q) const {
        std::vector<RigidBody*> out;
        for (auto* s : spheres) if (aabbOfSphere(*s).overlaps(q)) out.push_back(s->body);
        for (auto* b : boxes)   if (aabbOfBox(*b).overlaps(q))    out.push_back(b->body);
        return out;
    }
};

} // namespace phys
