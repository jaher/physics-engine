// The full joint zoo (hinge/revolute with limits + motor, slider/prismatic with
// limits + motor, fixed, distance, universal, gear) solved by iterative
// inverse-mass-weighted position projection (PBD style) with velocity-level
// motors — the constraint set exposed by ODE / MuJoCo / PhysX / Bullet / Unity.
#pragma once
#include "body.h"
#include <vector>

namespace phys {

// --- shared correction helpers ------------------------------------------------
// Move two bodies so worldA coincides with worldB, weighted by inverse mass and
// inverse inertia along the error direction (angular via anchor lever arms).
inline void projectPointCoincide(RigidBody* b0, RigidBody* b1, const Vector3& worldA, const Vector3& worldB, real slop = 0) {
    Vector3 delta = worldB - worldA; real len = delta.magnitude();
    if (len <= slop) return;
    Vector3 n = delta * (((real)1) / len);
    real mag = len - slop;
    RigidBody* bodies[2] = {b0, b1}; Vector3 anchors[2] = {worldA, worldB};
    real linI[2] = {0, 0}, angI[2] = {0, 0}, total = 0;
    for (int i = 0; i < 2; i++) { RigidBody* b = bodies[i]; if (!b || b->getInverseMass() == 0) continue;
        Vector3 rel = anchors[i] - b->getPosition();
        Matrix3 iit; b->getInverseInertiaTensorWorld(&iit);
        Vector3 ang = ((iit.transform(rel % n)) % rel);
        linI[i] = b->getInverseMass(); angI[i] = ang * n;
        total += linI[i] + angI[i]; }
    if (total <= 0) return;
    for (int i = 0; i < 2; i++) { RigidBody* b = bodies[i]; if (!b || b->getInverseMass() == 0) continue;
        real sign = i == 0 ? 1 : -1;
        real linMove = sign * mag * (linI[i] / total), angMove = sign * mag * (angI[i] / total);
        Vector3 pos = b->getPosition(); pos.addScaledVector(n, linMove); b->setPosition(pos);
        if (angI[i] > 1e-12) { Vector3 rel = anchors[i] - b->getPosition();
            Matrix3 iit; b->getInverseInertiaTensorWorld(&iit);
            Vector3 rot = iit.transform(rel % n) * (angMove / angI[i]);
            Quaternion q = b->getOrientation(); q.addScaledVector(rot, 1); b->setOrientation(q); }
        b->calculateDerivedData(); b->setAwake(true); }
}
// Apply a small relative rotation `rotvec` (axis*angle) between two bodies,
// weighted by inverse inertia about that axis.
inline void projectRelativeRotation(RigidBody* b0, RigidBody* b1, const Vector3& rotvec) {
    real ang = rotvec.magnitude(); if (ang < 1e-9) return;
    Vector3 axis = rotvec * (((real)1) / ang);
    real i0 = 0, i1 = 0;
    Matrix3 iit;
    if (b0 && b0->getInverseMass() > 0) { b0->getInverseInertiaTensorWorld(&iit); i0 = (iit.transform(axis)) * axis; }
    if (b1 && b1->getInverseMass() > 0) { b1->getInverseInertiaTensorWorld(&iit); i1 = (iit.transform(axis)) * axis; }
    real total = i0 + i1; if (total <= 0) return;
    if (i0 > 0) { Quaternion q = b0->getOrientation(); q.addScaledVector(axis * (ang * i0 / total), 1);
        b0->setOrientation(q); b0->calculateDerivedData(); b0->setAwake(true); }
    if (i1 > 0) { Quaternion q = b1->getOrientation(); q.addScaledVector(axis * (-ang * i1 / total), 1);
        b1->setOrientation(q); b1->calculateDerivedData(); b1->setAwake(true); }
}
// Velocity-level angular impulse about an axis: drives relative angular velocity
// (w0-w1)·axis toward `target`, with impulse magnitude capped by `maxImpulse`.
inline void angularVelocityMotor(RigidBody* b0, RigidBody* b1, const Vector3& axis, real target, real maxImpulse) {
    real rel = 0;
    if (b0) rel += b0->getRotation() * axis;
    if (b1) rel -= b1->getRotation() * axis;
    real i0 = 0, i1 = 0; Matrix3 iit;
    if (b0 && b0->getInverseMass() > 0) { b0->getInverseInertiaTensorWorld(&iit); i0 = (iit.transform(axis)) * axis; }
    if (b1 && b1->getInverseMass() > 0) { b1->getInverseInertiaTensorWorld(&iit); i1 = (iit.transform(axis)) * axis; }
    real total = i0 + i1; if (total <= 0) return;
    real lambda = (target - rel) / total;
    if (lambda > maxImpulse) lambda = maxImpulse; if (lambda < -maxImpulse) lambda = -maxImpulse;
    if (b0 && b0->getInverseMass() > 0) { b0->getInverseInertiaTensorWorld(&iit); b0->addRotation(iit.transform(axis * lambda)); b0->setAwake(true); }
    if (b1 && b1->getInverseMass() > 0) { b1->getInverseInertiaTensorWorld(&iit); b1->addRotation(iit.transform(axis * -lambda)); b1->setAwake(true); }
}

// --- joints --------------------------------------------------------------------
class Joint2 {
public:
    RigidBody* body[2] = {nullptr, nullptr};              // body[1] == nullptr → world
    virtual ~Joint2() {}
    virtual void project(real dt) = 0;                    // one solver iteration
};

class BallSocketJoint : public Joint2 {
public:
    Vector3 local[2];                                     // anchor in each body's space
    void set(RigidBody* a, const Vector3& la, RigidBody* b, const Vector3& lb) {
        body[0] = a; body[1] = b; local[0] = la; local[1] = lb;
    }
    Vector3 worldAnchor(int i) const { return body[i] ? body[i]->getPointInWorldSpace(local[i]) : local[i]; }
    void project(real) override { projectPointCoincide(body[0], body[1], worldAnchor(0), worldAnchor(1)); }
};

class DistanceJoint : public Joint2 {
public:
    Vector3 local[2]; real length = 1;
    void set(RigidBody* a, const Vector3& la, RigidBody* b, const Vector3& lb, real len) {
        body[0] = a; body[1] = b; local[0] = la; local[1] = lb; length = len;
    }
    void project(real) override {
        Vector3 wa = body[0] ? body[0]->getPointInWorldSpace(local[0]) : local[0];
        Vector3 wb = body[1] ? body[1]->getPointInWorldSpace(local[1]) : local[1];
        Vector3 d = wb - wa; real len = d.magnitude(); if (len < 1e-9) return;
        Vector3 target = wb - d * (length / len);         // point at correct distance
        projectPointCoincide(body[0], body[1], wa, target);
    }
};

class HingeJoint : public Joint2 {
public:
    Vector3 localAnchor[2], localAxis[2];
    Vector3 localRef0, localRef1;                          // perpendicular reference dirs for the angle
    bool hasLimits = false; real minAngle = 0, maxAngle = 0;
    bool hasMotor = false; real motorTargetVel = 0, motorMaxImpulse = (real)1e9;
    void set(RigidBody* a, const Vector3& anchorA, const Vector3& axisA,
             RigidBody* b, const Vector3& anchorB, const Vector3& axisB) {
        body[0] = a; body[1] = b;
        localAnchor[0] = anchorA; localAnchor[1] = anchorB;
        localAxis[0] = axisA.unit(); localAxis[1] = axisB.unit();
        // build reference vectors perpendicular to the axis for angle measurement
        Vector3 t = real_abs(localAxis[0].x) < (real)0.9 ? Vector3(1, 0, 0) : Vector3(0, 1, 0);
        localRef0 = (localAxis[0] % t).unit();
        Vector3 worldRef = a ? a->getDirectionInWorldSpace(localRef0) : localRef0;
        localRef1 = b ? b->getDirectionInLocalSpace(worldRef) : worldRef;
    }
    Vector3 worldAxis(int i) const { return body[i] ? body[i]->getDirectionInWorldSpace(localAxis[i]) : localAxis[i]; }
    real angle() const {
        Vector3 a0 = worldAxis(0);
        Vector3 r0 = body[0] ? body[0]->getDirectionInWorldSpace(localRef0) : localRef0;
        Vector3 r1 = body[1] ? body[1]->getDirectionInWorldSpace(localRef1) : localRef1;
        Vector3 r1p = (r1 - a0 * (r1 * a0)).unit();
        real c = r0 * r1p, s = (r0 % r1p) * a0;
        return std::atan2((double)s, (double)c);
    }
    void project(real dt) override {
        Vector3 wa = body[0] ? body[0]->getPointInWorldSpace(localAnchor[0]) : localAnchor[0];
        Vector3 wb = body[1] ? body[1]->getPointInWorldSpace(localAnchor[1]) : localAnchor[1];
        projectPointCoincide(body[0], body[1], wa, wb);
        Vector3 a0 = worldAxis(0), a1 = worldAxis(1);
        projectRelativeRotation(body[0], body[1], a0 % a1);   // align the axes
        if (hasLimits) { real ang = angle();
            if (ang < minAngle) projectRelativeRotation(body[0], body[1], worldAxis(0) * (ang - minAngle));
            else if (ang > maxAngle) projectRelativeRotation(body[0], body[1], worldAxis(0) * (ang - maxAngle)); }
        if (hasMotor) angularVelocityMotor(body[0], body[1], worldAxis(0), motorTargetVel, motorMaxImpulse * dt);
    }
};

class SliderJoint : public Joint2 {
public:
    Vector3 localAnchor[2], localAxis[2];
    Quaternion refRel;                                     // locked relative orientation
    bool hasLimits = false; real minT = 0, maxT = 0;
    bool hasMotor = false; real motorTargetVel = 0, motorMaxImpulse = (real)1e9;
    void set(RigidBody* a, const Vector3& anchorA, const Vector3& axisA,
             RigidBody* b, const Vector3& anchorB, const Vector3& axisB) {
        body[0] = a; body[1] = b;
        localAnchor[0] = anchorA; localAnchor[1] = anchorB;
        localAxis[0] = axisA.unit(); localAxis[1] = axisB.unit();
    }
    Vector3 worldAxis() const { return body[0] ? body[0]->getDirectionInWorldSpace(localAxis[0]) : localAxis[0]; }
    real translation() const {                              // +ve when body[1] sits along +axis from body[0]
        Vector3 wa = body[0] ? body[0]->getPointInWorldSpace(localAnchor[0]) : localAnchor[0];
        Vector3 wb = body[1] ? body[1]->getPointInWorldSpace(localAnchor[1]) : localAnchor[1];
        return (wb - wa) * worldAxis();
    }
    void project(real dt) override {
        Vector3 a0 = worldAxis();
        Vector3 a1 = body[1] ? body[1]->getDirectionInWorldSpace(localAxis[1]) : localAxis[1];
        projectRelativeRotation(body[0], body[1], a0 % a1);   // keep axes parallel (no relative twist about ⊥)
        Vector3 wa = body[0] ? body[0]->getPointInWorldSpace(localAnchor[0]) : localAnchor[0];
        Vector3 wb = body[1] ? body[1]->getPointInWorldSpace(localAnchor[1]) : localAnchor[1];
        Vector3 d = wb - wa;
        // kill off-axis displacement: body1's anchor coincides with its rail projection
        projectPointCoincide(body[0], body[1], wa + a0 * (d * a0), wb);
        if (hasLimits) { real t = translation();
            if (t < minT) projectPointCoincide(body[0], body[1], wb + a0 * (minT - t), wb);
            else if (t > maxT) projectPointCoincide(body[0], body[1], wb - a0 * (t - maxT), wb); }
        if (hasMotor) {                                     // drives d(translation)/dt = (v1-v0)·axis → target
            RigidBody* b0 = body[0]; RigidBody* b1 = body[1];
            real rel = 0;
            if (b1) rel += b1->getVelocity() * a0; if (b0) rel -= b0->getVelocity() * a0;
            real m0 = b0 ? b0->getInverseMass() : 0, m1 = b1 ? b1->getInverseMass() : 0;
            real total = m0 + m1; if (total <= 0) return;
            real lam = (motorTargetVel - rel) / total;
            real cap = motorMaxImpulse * dt; if (lam > cap) lam = cap; if (lam < -cap) lam = -cap;
            if (m1 > 0) b1->addVelocity(a0 * (lam * m1));
            if (m0 > 0) b0->addVelocity(a0 * (-lam * m0));
        }
    }
};

class FixedJoint : public Joint2 {
public:
    Vector3 local[2]; Quaternion refRel;                   // b1-relative orientation at attach time
    void set(RigidBody* a, const Vector3& la, RigidBody* b, const Vector3& lb) {
        body[0] = a; body[1] = b; local[0] = la; local[1] = lb;
        Quaternion qa = a ? a->getOrientation() : Quaternion();
        Quaternion qb = b ? b->getOrientation() : Quaternion();
        Quaternion qbInv(qb.r, -qb.i, -qb.j, -qb.k);
        refRel = qbInv; refRel *= qa;                       // qb⁻¹ · qa
    }
    void project(real) override {
        Vector3 wa = body[0] ? body[0]->getPointInWorldSpace(local[0]) : local[0];
        Vector3 wb = body[1] ? body[1]->getPointInWorldSpace(local[1]) : local[1];
        projectPointCoincide(body[0], body[1], wa, wb);
        Quaternion qa = body[0] ? body[0]->getOrientation() : Quaternion();
        Quaternion qb = body[1] ? body[1]->getOrientation() : Quaternion();
        Quaternion target = qb; target *= refRel;           // desired qa
        // rotation error qe = target⁻¹ · qa  → small-angle rotvec
        Quaternion ti(target.r, -target.i, -target.j, -target.k);
        Quaternion qe = ti; qe *= qa;
        if (qe.r < 0) { qe.r = -qe.r; qe.i = -qe.i; qe.j = -qe.j; qe.k = -qe.k; }
        Vector3 rotvec(qe.i * 2, qe.j * 2, qe.k * 2);       // ≈ axis*angle for small errors
        Vector3 world = body[0] ? body[0]->getDirectionInWorldSpace(rotvec) : rotvec;
        projectRelativeRotation(body[0], body[1], world * -1);
    }
};

// Universal joint: point coincidence + the two shaft axes stay perpendicular.
class UniversalJoint : public Joint2 {
public:
    Vector3 localAnchor[2], localAxis[2];
    void set(RigidBody* a, const Vector3& anchorA, const Vector3& axisA,
             RigidBody* b, const Vector3& anchorB, const Vector3& axisB) {
        body[0] = a; body[1] = b; localAnchor[0] = anchorA; localAnchor[1] = anchorB;
        localAxis[0] = axisA.unit(); localAxis[1] = axisB.unit();
    }
    void project(real) override {
        Vector3 wa = body[0] ? body[0]->getPointInWorldSpace(localAnchor[0]) : localAnchor[0];
        Vector3 wb = body[1] ? body[1]->getPointInWorldSpace(localAnchor[1]) : localAnchor[1];
        projectPointCoincide(body[0], body[1], wa, wb);
        Vector3 a0 = body[0] ? body[0]->getDirectionInWorldSpace(localAxis[0]) : localAxis[0];
        Vector3 a1 = body[1] ? body[1]->getDirectionInWorldSpace(localAxis[1]) : localAxis[1];
        real dot = a0 * a1;                                 // want 0 (perpendicular)
        Vector3 corrAxis = a0 % a1; real cl = corrAxis.magnitude(); if (cl < 1e-9) return;
        projectRelativeRotation(body[0], body[1], corrAxis * (dot / cl));
    }
};

// Gear: couples two hinge axes so w0·a0 = -ratio · (w1·a1)  (velocity level).
class GearJoint : public Joint2 {
public:
    Vector3 localAxis[2]; real ratio = 1;
    void set(RigidBody* a, const Vector3& axisA, RigidBody* b, const Vector3& axisB, real ratio_) {
        body[0] = a; body[1] = b; localAxis[0] = axisA.unit(); localAxis[1] = axisB.unit(); ratio = ratio_;
    }
    void project(real) override {
        Vector3 a0 = body[0] ? body[0]->getDirectionInWorldSpace(localAxis[0]) : localAxis[0];
        Vector3 a1 = body[1] ? body[1]->getDirectionInWorldSpace(localAxis[1]) : localAxis[1];
        real w0 = body[0] ? body[0]->getRotation() * a0 : 0;
        real w1 = body[1] ? body[1]->getRotation() * a1 : 0;
        real err = w0 + ratio * w1;                        // want 0
        real i0 = 0, i1 = 0; Matrix3 iit;
        if (body[0] && body[0]->getInverseMass() > 0) { body[0]->getInverseInertiaTensorWorld(&iit); i0 = (iit.transform(a0)) * a0; }
        if (body[1] && body[1]->getInverseMass() > 0) { body[1]->getInverseInertiaTensorWorld(&iit); i1 = (iit.transform(a1)) * a1; }
        real total = i0 + ratio * ratio * i1; if (total <= 0) return;
        real lam = -err / total;
        if (i0 > 0) { body[0]->getInverseInertiaTensorWorld(&iit); body[0]->addRotation(iit.transform(a0 * lam)); }
        if (i1 > 0) { body[1]->getInverseInertiaTensorWorld(&iit); body[1]->addRotation(iit.transform(a1 * (lam * ratio))); }
    }
};

// PD position servo on a hinge (MuJoCo-style actuator): drives angle → target.
class HingeServo {
public:
    HingeJoint* hinge = nullptr; real target = 0, kp = 8, kd = 1, maxVel = 20;
    void update(real dt) {
        if (!hinge) return;
        real err = target - hinge->angle();
        real vel = kp * err;                               // kd folded into the motor's implicit damping
        if (vel > maxVel) vel = maxVel; if (vel < -maxVel) vel = -maxVel;
        hinge->hasMotor = true; hinge->motorTargetVel = vel;
        (void)dt; (void)kd;
    }
};

class JointSolver {
public:
    std::vector<Joint2*> joints;
    int iterations = 8;
    void solve(real dt) { for (int it = 0; it < iterations; it++) for (auto* j : joints) j->project(dt); }
};

} // namespace phys
