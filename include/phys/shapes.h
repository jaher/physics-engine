// Capsule collision primitive + narrow-phase detectors (capsule-plane,
// capsule-sphere, capsule-capsule, capsule-box). Closes the shape gap with
// MuJoCo/PhysX/Bullet, whose workhorse primitive is the capsule.
#pragma once
#include "collide_fine.h"

namespace phys {

// Capsule along its local Y axis: segment of half-length `halfHeight`, radius r.
class CollisionCapsule : public CollisionPrimitive {
public:
    real radius = (real)0.5, halfHeight = (real)0.5;
    void endpoints(Vector3& a, Vector3& b) const {
        Vector3 c = getAxis(3), ax = getAxis(1);
        a = c + ax * halfHeight; b = c - ax * halfHeight;
    }
};

inline Vector3 closestOnSegment(const Vector3& a, const Vector3& b, const Vector3& p) {
    Vector3 ab = b - a; real t = ((p - a) * ab) / (ab * ab + (real)1e-12);
    if (t < 0) t = 0; if (t > 1) t = 1;
    return a + ab * t;
}
// closest points between two segments (standard clamped solve)
inline void closestSegSeg(const Vector3& p1, const Vector3& q1, const Vector3& p2, const Vector3& q2,
                          Vector3& c1, Vector3& c2) {
    Vector3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    real a = d1 * d1, e = d2 * d2, f = d2 * r, s, t;
    if (a <= 1e-12 && e <= 1e-12) { c1 = p1; c2 = p2; return; }
    if (a <= 1e-12) { s = 0; t = f / e; t = t < 0 ? 0 : (t > 1 ? 1 : t); }
    else { real c = d1 * r;
        if (e <= 1e-12) { t = 0; s = -c / a; s = s < 0 ? 0 : (s > 1 ? 1 : s); }
        else { real b = d1 * d2, den = a * e - b * b;
            s = den != 0 ? (b * f - c * e) / den : 0; s = s < 0 ? 0 : (s > 1 ? 1 : s);
            t = (b * s + f) / e;
            if (t < 0) { t = 0; s = -c / a; s = s < 0 ? 0 : (s > 1 ? 1 : s); }
            else if (t > 1) { t = 1; s = (b - c) / a; s = s < 0 ? 0 : (s > 1 ? 1 : s); } } }
    c1 = p1 + d1 * s; c2 = p2 + d2 * t;
}

struct CapsuleDetector {
    static unsigned capsuleAndHalfSpace(const CollisionCapsule& cap, const CollisionPlane& plane, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 e[2]; cap.endpoints(e[0], e[1]);
        unsigned used = 0;
        for (int i = 0; i < 2 && data->hasMoreContacts(); i++) {
            real dist = plane.direction * e[i] - cap.radius - plane.offset;
            if (dist >= 0) continue;
            Contact* c = data->contacts;
            c->contactNormal = plane.direction;
            c->penetration = -dist;
            c->contactPoint = e[i] - plane.direction * cap.radius;
            c->setBodyData(cap.body, nullptr, data->friction, data->restitution);
            data->addContacts(1); used++;
        }
        return used;
    }
    static unsigned capsuleAndSphere(const CollisionCapsule& cap, const CollisionSphere& sph, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 a, b; cap.endpoints(a, b);
        Vector3 sc = sph.getAxis(3);
        Vector3 p = closestOnSegment(a, b, sc);
        Vector3 mid = p - sc; real dist = mid.magnitude();
        if (dist <= 0 || dist >= cap.radius + sph.radius) return 0;
        Contact* c = data->contacts;
        c->contactNormal = mid * (((real)1) / dist);            // from sphere toward capsule
        c->penetration = cap.radius + sph.radius - dist;
        c->contactPoint = sc + mid * (sph.radius / (cap.radius + sph.radius));
        c->setBodyData(cap.body, sph.body, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    static unsigned capsuleAndCapsule(const CollisionCapsule& one, const CollisionCapsule& two, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 a1, b1, a2, b2; one.endpoints(a1, b1); two.endpoints(a2, b2);
        Vector3 c1, c2; closestSegSeg(a1, b1, a2, b2, c1, c2);
        Vector3 mid = c1 - c2; real dist = mid.magnitude();
        real rsum = one.radius + two.radius;
        if (dist <= 0 || dist >= rsum) return 0;
        Contact* c = data->contacts;
        c->contactNormal = mid * (((real)1) / dist);
        c->penetration = rsum - dist;
        c->contactPoint = c2 + mid * (two.radius / rsum);
        c->setBodyData(one.body, two.body, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    // capsule vs box: treat the capsule's closest segment point to the box as a sphere
    static unsigned capsuleAndBox(const CollisionCapsule& cap, const CollisionBox& box, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 a, b; cap.endpoints(a, b);
        // find segment point closest to box centre in box space, then clamp to box
        unsigned used = 0;
        for (int pass = 0; pass < 2 && data->hasMoreContacts(); pass++) {
            Vector3 probe = pass == 0 ? a : b;
            Vector3 rel = box.getTransform().transformInverse(probe);
            Vector3 clamped = rel;
            for (int i = 0; i < 3; i++) { if (clamped[i] > box.halfSize[i]) clamped[i] = box.halfSize[i];
                if (clamped[i] < -box.halfSize[i]) clamped[i] = -box.halfSize[i]; }
            Vector3 closestWorld = box.getTransform().transform(clamped);
            Vector3 seg = closestOnSegment(a, b, closestWorld);
            Vector3 d = seg - closestWorld; real dist = d.magnitude();
            if (dist <= 0 || dist >= cap.radius) continue;
            Contact* c = data->contacts;
            c->contactNormal = d * (((real)1) / dist);
            c->penetration = cap.radius - dist;
            c->contactPoint = closestWorld;
            c->setBodyData(cap.body, box.body, data->friction, data->restitution);
            data->addContacts(1); used++;
        }
        return used;
    }
};

} // namespace phys
