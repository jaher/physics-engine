// Grab-bag of engine features from the survey: continuous collision detection
// (swept spheres), contact events, trigger volumes, kinematic bodies, and an RK4
// integrator (Brax/Chrono-style) for force-field particles.
#pragma once
#include "raycast.h"
#include <functional>
#include <set>

namespace phys {

// ---- continuous collision detection (CCD) ------------------------------------
// Earliest time t in [0,1] a sphere sweeping p0→p1 hits a plane / static sphere.
inline bool sweepSpherePlane(const Vector3& p0, const Vector3& p1, real r,
                             const Vector3& n, real offset, real& t, Vector3& hitN) {
    real d0 = n * p0 - offset, d1 = n * p1 - offset;
    if (d0 >= r && d1 >= r) return false;                  // stays clear
    if (d0 < r) { t = 0; hitN = n; return true; }          // already touching
    t = (d0 - r) / (d0 - d1); hitN = n; return true;
}
inline bool sweepSphereSphere(const Vector3& p0, const Vector3& p1, real r,
                              const Vector3& c, real cr, real& t, Vector3& hitN) {
    Vector3 d = p1 - p0, m = p0 - c; real R = r + cr;
    real a = d * d; if (a < 1e-12) return false;
    real b = 2 * (m * d), cc = m * m - R * R;
    if (cc < 0) { t = 0; hitN = (p0 - c).unit(); return true; }
    real disc = b * b - 4 * a * cc; if (disc < 0) return false;
    real tt = (-b - real_sqrt(disc)) / (2 * a);
    if (tt < 0 || tt > 1) return false;
    t = tt; hitN = ((p0 + d * tt) - c).unit(); return true;
}
// Integrate a fast body's position with tunnel-proof sweeps against a set of
// planes and static spheres; on hit, stop at the surface and reflect velocity.
struct CCDWorld {
    struct SPlane { Vector3 n; real offset; };
    struct SSphere { Vector3 c; real r; };
    std::vector<SPlane> planes; std::vector<SSphere> spheres;
    // returns true if a hit occurred this step
    bool integrate(RigidBody& b, real radius, real dt, real restitution) {
        Vector3 p0 = b.getPosition();
        b.integrate(dt);
        Vector3 p1 = b.getPosition();
        real bestT = 2; Vector3 bestN;
        real t; Vector3 n;
        for (auto& pl : planes) if (sweepSpherePlane(p0, p1, radius, pl.n, pl.offset, t, n) && t < bestT) { bestT = t; bestN = n; }
        for (auto& sp : spheres) if (sweepSphereSphere(p0, p1, radius, sp.c, sp.r, t, n) && t < bestT) { bestT = t; bestN = n; }
        if (bestT > 1) return false;
        Vector3 hit = p0 + (p1 - p0) * bestT;
        b.setPosition(hit + bestN * (real)1e-4);
        Vector3 v = b.getVelocity();
        real vn = v * bestN;
        if (vn < 0) b.setVelocity(v - bestN * ((1 + restitution) * vn));
        b.calculateDerivedData();
        return true;
    }
};

// ---- contact events (Unity/PhysX OnCollisionEnter/Exit) -----------------------
class ContactEvents {
public:
    using Pair = std::pair<RigidBody*, RigidBody*>;
    std::function<void(RigidBody*, RigidBody*)> onBegin, onEnd;
    void frame(const Contact* contacts, unsigned count) {
        std::set<Pair> now;
        for (unsigned i = 0; i < count; i++) {
            Pair p(contacts[i].body[0], contacts[i].body[1]);
            if (p.first > p.second) std::swap(p.first, p.second);
            now.insert(p);
        }
        if (onBegin) for (auto& p : now) if (!prev.count(p)) onBegin(p.first, p.second);
        if (onEnd) for (auto& p : prev) if (!now.count(p)) onEnd(p.first, p.second);
        prev = now;
    }
private:
    std::set<Pair> prev;
};

// ---- trigger volumes (no collision response, just enter/exit callbacks) -------
class TriggerVolume {
public:
    Vector3 centre; Vector3 halfSize; bool isSphere = false; real radius = 1;
    std::function<void(RigidBody*)> onEnter, onExit;
    void frame(const std::vector<RigidBody*>& bodies) {
        std::set<RigidBody*> now;
        for (auto* b : bodies) {
            Vector3 p = b->getPosition() - centre;
            bool inside = isSphere ? (p.magnitude() <= radius)
                : (real_abs(p.x) <= halfSize.x && real_abs(p.y) <= halfSize.y && real_abs(p.z) <= halfSize.z);
            if (inside) now.insert(b);
        }
        if (onEnter) for (auto* b : now) if (!prev.count(b)) onEnter(b);
        if (onExit) for (auto* b : prev) if (!now.count(b)) onExit(b);
        prev = now;
    }
private:
    std::set<RigidBody*> prev;
};

// ---- kinematic bodies (animated, infinite mass, push dynamics) ----------------
// Give the body infinite mass but a scripted velocity: contacts then see the
// motion and push dynamic bodies correctly.
inline void makeKinematic(RigidBody& b) {
    b.setInverseMass(0); b.setInverseInertiaTensor(Matrix3());
    b.setCanSleep(false); b.setAwake(true); b.setAcceleration(0, 0, 0);
}
inline void moveKinematic(RigidBody& b, const Vector3& velocity, real dt) {
    b.setVelocity(velocity);
    b.setPosition(b.getPosition() + velocity * dt);
    b.calculateDerivedData();
}

// ---- RK4 integrator (Brax/Chrono offer higher-order integration) ---------------
// Integrates x'' = accel(x, v, t) one step with classical Runge-Kutta 4.
inline void integrateRK4(Vector3& x, Vector3& v, real t, real dt,
                         const std::function<Vector3(const Vector3&, const Vector3&, real)>& accel) {
    auto f = [&](const Vector3& px, const Vector3& pv, real pt, Vector3& dx, Vector3& dv) {
        dx = pv; dv = accel(px, pv, pt);
    };
    Vector3 k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;
    f(x, v, t, k1x, k1v);
    f(x + k1x * (dt / 2), v + k1v * (dt / 2), t + dt / 2, k2x, k2v);
    f(x + k2x * (dt / 2), v + k2v * (dt / 2), t + dt / 2, k3x, k3v);
    f(x + k3x * dt, v + k3v * dt, t + dt, k4x, k4v);
    x += (k1x + (k2x + k3x) * 2 + k4x) * (dt / 6);
    v += (k1v + (k2v + k3v) * 2 + k4v) * (dt / 6);
}

} // namespace phys
