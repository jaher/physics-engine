// The robotics layer this engine was missing: sensors (IMU, Lidar, contact-force),
// recursive Newton–Euler inverse dynamics for a serial revolute arm, a
// Jacobian-transpose inverse-kinematics solver, and tendon/muscle actuators.
// Built only on core.h (Vector3/Matrix3/Quaternion), so it needs no rigid-body
// state — the sensors are fed measured quantities, the arm is its own minimal
// serial-chain struct. References: Craig, *Introduction to Robotics*, ch.6 (RNE);
// Buss, *Introduction to Inverse Kinematics* (Jacobian transpose); Hill (muscle).
#pragma once
#include "core.h"
#include <vector>

namespace phys {

// ============================================================ (1) sensors =====

// IMU: a strapdown accelerometer + gyro rigidly mounted on a body. The
// accelerometer measures *specific force* (proper acceleration) a − g expressed
// in the body frame, so a body at rest reads g pointing "up" (the reaction that
// holds it against gravity); the gyro reports the body-frame angular rate.
struct IMU {
    Vector3 gravity = Vector3(0, -9.81, 0);                 // world gravitational acceleration
    struct Reading { Vector3 accel, gyro; };
    // linAccWorld/angVelWorld: the body's true linear accel & angular velocity in
    // world; orientation: body→world. Both outputs are in the (body/sensor) frame.
    Reading measure(const Vector3& linAccWorld, const Vector3& angVelWorld,
                    const Quaternion& orientation) const {
        Matrix3 R; R.setOrientation(orientation);           // body→world
        Reading r;
        r.accel = R.transformTranspose(linAccWorld - gravity);   // world→body: proper accel
        r.gyro  = R.transformTranspose(angVelWorld);             // world→body: angular rate
        return r;
    }
};

// Lidar: a planar fan of N rays swept over `fov`, each range resolved by a
// user-supplied raycast callable `real(Vector3 origin, Vector3 dir)` (returns the
// hit distance). The fan lies in the plane spanned by `forward` and `up × forward`.
struct Lidar {
    int numRays = 16;
    real fov = real_pi;                                     // total angular span
    real minRange = 0, maxRange = 100;
    Vector3 origin;
    Vector3 forward = Vector3(1, 0, 0), up = Vector3(0, 1, 0);

    template <class Raycast>
    std::vector<real> scan(Raycast raycast) const {
        std::vector<real> ranges((size_t)(numRays > 0 ? numRays : 0));
        Vector3 f = forward.unit();
        Vector3 left = (up.unit() % f).unit();              // in-plane, ⟂ to forward
        for (int i = 0; i < numRays; i++) {
            real a = (numRays == 1) ? (real)0
                                    : -fov * (real)0.5 + fov * ((real)i / (numRays - 1));
            Vector3 dir = (f * real_cos(a) + left * real_sin(a)).unit();
            real r = raycast(origin, dir);
            if (r < minRange) r = minRange;
            if (r > maxRange) r = maxRange;
            ranges[(size_t)i] = r;
        }
        return ranges;
    }
};

// ContactForceSensor: integrate the contact impulses landing on a body over a
// step and divide by dt to recover the average contact force (force = ∫J / dt).
struct ContactForceSensor {
    Vector3 impulseAccum;
    void addImpulse(const Vector3& j) { impulseAccum += j; }
    void clear() { impulseAccum.clear(); }
    Vector3 force(real dt) const { return dt > 0 ? impulseAccum * (((real)1) / dt) : Vector3(); }
};

// ============================================ (2) serial chain + RNE / IK =====

// One revolute link of a serial arm, all vectors in the link's own frame:
//   axis    — the joint rotation axis
//   toChild — offset from this joint to the next joint (the link's length)
//   com     — offset from this joint to the link's centre of mass
//   inertia — inertia tensor about the COM (zero ⇒ point mass)
struct ChainLink {
    Vector3 axis    = Vector3(0, 0, 1);
    Vector3 toChild = Vector3(1, 0, 0);
    Vector3 com     = Vector3(0.5, 0, 0);
    real    mass    = 1;
    Matrix3 inertia;                                        // defaults to all-zero
};

// Forward-kinematics snapshot: world-frame axes, joint positions (n+1: joints
// plus the tip), COM positions, joint→joint vectors and world inertia tensors.
struct FKResult {
    std::vector<Vector3> jointPos;                          // size n+1
    std::vector<Vector3> axis, comPos, dvec;                // size n
    std::vector<Matrix3> Iworld;                            // size n
    Vector3 endEffector;
};

struct SerialChain {
    std::vector<ChainLink> links;
    Vector3 gravity = Vector3(0, -9.81, 0);
    Vector3 baseOrigin;
    Quaternion baseOrient;

    static Quaternion axisAngle(const Vector3& axis, real angle) {
        Vector3 a = axis.unit();
        real h = angle * (real)0.5, s = real_sin(h);
        return Quaternion(real_cos(h), a.x * s, a.y * s, a.z * s);
    }

    // Walk the chain, composing each joint's rotation about its axis, to place
    // every link frame in the world.
    FKResult fk(const std::vector<real>& q) const {
        FKResult r;
        int n = (int)links.size();
        Matrix3 Rparent; Rparent.setOrientation(baseOrient);
        Vector3 Pjoint = baseOrigin;
        r.jointPos.push_back(Pjoint);
        for (int i = 0; i < n; i++) {
            const ChainLink& L = links[i];
            Vector3 zw = (Rparent * L.axis).unit();          // world joint axis (invariant under its own spin)
            Matrix3 Rj; Rj.setOrientation(axisAngle(L.axis, q[(size_t)i]));
            Matrix3 Ri = Rparent * Rj;                        // link i world orientation
            Vector3 comW  = Pjoint + Ri * L.com;
            Vector3 d     = Ri * L.toChild;
            Vector3 Pnext = Pjoint + d;
            r.axis.push_back(zw);
            r.comPos.push_back(comW);
            r.dvec.push_back(d);
            r.Iworld.push_back(Ri * L.inertia * Ri.transpose());
            r.jointPos.push_back(Pnext);
            Rparent = Ri; Pjoint = Pnext;
        }
        r.endEffector = Pjoint;
        return r;
    }

    Vector3 endEffector(const std::vector<real>& q) const { return fk(q).endEffector; }

    // Recursive Newton–Euler inverse dynamics: joint torques that realise the
    // desired (q, qd, qdd). Gravity enters as a fictitious base acceleration −g.
    // An outward pass propagates velocity/acceleration; an inward pass gathers
    // forces/moments and projects each onto its joint axis.
    std::vector<real> inverseDynamics(const std::vector<real>& q,
                                      const std::vector<real>& qd,
                                      const std::vector<real>& qdd) const {
        int n = (int)links.size();
        FKResult k = fk(q);
        std::vector<Vector3> omega(n), alpha(n), aJoint(n + 1), aCom(n), F(n), N(n), f(n + 1), nn(n + 1);
        aJoint[0] = -gravity;                                // base origin "accelerates" at −g
        for (int i = 0; i < n; i++) {                        // outward: kinematics
            Vector3 z = k.axis[i];
            Vector3 wPrev = (i == 0) ? Vector3() : omega[i - 1];
            Vector3 aPrev = (i == 0) ? Vector3() : alpha[i - 1];
            omega[i] = wPrev + z * qd[(size_t)i];
            alpha[i] = aPrev + z * qdd[(size_t)i] + (wPrev % (z * qd[(size_t)i]));
            Vector3 rc = k.comPos[i] - k.jointPos[i];
            aCom[i]      = aJoint[i] + (alpha[i] % rc)        + (omega[i] % (omega[i] % rc));
            aJoint[i + 1]= aJoint[i] + (alpha[i] % k.dvec[i]) + (omega[i] % (omega[i] % k.dvec[i]));
        }
        std::vector<real> tau((size_t)n);
        for (int i = n - 1; i >= 0; i--) {                   // inward: forces & moments
            F[i] = aCom[i] * links[(size_t)i].mass;
            N[i] = k.Iworld[i] * alpha[i] + (omega[i] % (k.Iworld[i] * omega[i]));
            f[i] = F[i] + f[i + 1];
            Vector3 rc = k.comPos[i] - k.jointPos[i];
            nn[i] = N[i] + nn[i + 1] + (rc % F[i]) + (k.dvec[i] % f[i + 1]);
            tau[(size_t)i] = nn[i] * k.axis[i];               // project the joint moment onto its axis
        }
        return tau;
    }
};

// Jacobian-transpose IK: nudge q so the end-effector approaches `target`. Each
// column of the position Jacobian is zᵢ × (p_ee − pᵢ); the step is a scaled
// Jᵀe, with Buss's optimal scalar (× a safety factor) so the error decreases.
struct JacobianIK {
    const SerialChain* chain = nullptr;
    real safety  = (real)0.5;                               // shrink the optimal step to stay stable
    real maxStep = (real)0.3;                               // per-joint clamp (rad)

    real errorNorm(const std::vector<real>& q, const Vector3& target) const {
        return (target - chain->endEffector(q)).magnitude();
    }

    // One iteration; returns the error norm AFTER the update.
    real step(std::vector<real>& q, const Vector3& target) const {
        int n = (int)chain->links.size();
        FKResult k = chain->fk(q);
        Vector3 ee = k.endEffector, e = target - ee;
        std::vector<Vector3> J((size_t)n);
        std::vector<real>    JTe((size_t)n);
        Vector3 JJTe;
        for (int i = 0; i < n; i++) {
            J[i]   = k.axis[i] % (ee - k.jointPos[i]);       // position Jacobian column
            JTe[i] = J[i] * e;                                // (Jᵀe)_i
        }
        for (int i = 0; i < n; i++) JJTe += J[i] * JTe[i];
        real denom = JJTe * JJTe;
        real alpha = (denom > (real)1e-12) ? (e * JJTe) / denom : 0;
        alpha *= safety;
        for (int i = 0; i < n; i++) {
            real dq = alpha * JTe[i];
            if (dq >  maxStep) dq =  maxStep;
            if (dq < -maxStep) dq = -maxStep;
            q[(size_t)i] += dq;
        }
        return (target - chain->endEffector(q)).magnitude();
    }

    // Run up to `iters` steps; optionally record the error at each step.
    real solve(std::vector<real>& q, const Vector3& target, int iters,
               std::vector<real>* history = nullptr) const {
        real err = errorNorm(q, target);
        if (history) history->push_back(err);
        for (int it = 0; it < iters; it++) { err = step(q, target); if (history) history->push_back(err); }
        return err;
    }
};

// ============================================= (4) tendon / muscle actuators ==

// Tendon: a cable routed through world-space waypoints; its tension is a linear
// spring-damper on the routed length, one-sided (a cable can only pull, ≥ 0).
struct Tendon {
    std::vector<Vector3> route;
    real restLength = 0, stiffness = 1000, damping = 10;
    real length() const {
        real L = 0;
        for (size_t i = 1; i < route.size(); i++) L += (route[i] - route[i - 1]).magnitude();
        return L;
    }
    real tension(real len, real lenRate) const {
        real f = stiffness * (len - restLength) + damping * lenRate;
        return f > 0 ? f : 0;                                // cannot push
    }
    real tension() const { return tension(length(), 0); }
};

// HillMuscle: a Hill-type actuator — active force = activation·Fmax·f_L(L) with a
// parabolic force-length bell peaking at optLength, plus a one-sided passive
// spring beyond optLength and optional damping. Output tension is ≥ 0.
struct HillMuscle {
    real Fmax = 1000, optLength = (real)0.1, width = (real)0.5;   // width of the active bell (fraction of optLength)
    real activation = 1;                                          // neural drive 0..1
    real passiveStiffness = 0, damping = 0;
    real activeForceLength(real L) const {
        real x = (L - optLength) / (width * optLength);
        real f = 1 - x * x;
        return f > 0 ? f : 0;
    }
    real passiveForceLength(real L) const { real d = L - optLength; return d > 0 ? passiveStiffness * d : 0; }
    real tension(real L, real Ldot) const {
        real f = activation * Fmax * activeForceLength(L) + passiveForceLength(L) + damping * Ldot;
        return f > 0 ? f : 0;
    }
};

} // namespace phys
