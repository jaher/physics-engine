// A rigid body: position + orientation, linear + angular velocity, mass and an
// inertia tensor (Millington Part III ch.10). Integration advances all six
// degrees of freedom and maintains the derived world transform and world-space
// inverse inertia tensor. Includes the sleep (deactivation) system.
#pragma once
#include "core.h"
#include <cassert>

namespace phys {

class RigidBody {
public:
    static real sleepEpsilon;
protected:
    real inverseMass = 0;
    Matrix3 inverseInertiaTensor;              // body space
    real linearDamping = (real)0.99;
    real angularDamping = (real)0.99;
    Vector3 position;
    Quaternion orientation;
    Vector3 velocity;                          // linear
    Vector3 rotation;                          // angular velocity
    Matrix3 inverseInertiaTensorWorld;         // derived
    real motion = 0;                           // recency-weighted kinetic energy, for sleep
    bool isAwakeFlag = true;
    bool canSleepFlag = true;
    Matrix4 transformMatrix;                   // derived body→world
    Vector3 forceAccum, torqueAccum;
    Vector3 acceleration;                      // constant (e.g. gravity)
    Vector3 lastFrameAcceleration;             // derived, for contact resolution
public:
    RigidBody() { motion = sleepEpsilon * 2; }     // start above the sleep threshold

    // world inverse-inertia = R · I_body⁻¹ · Rᵀ
    static void transformInertiaTensor(Matrix3& iitWorld, const Matrix3& iitBody, const Matrix4& rot) {
        Matrix3 R(rot.data[0], rot.data[1], rot.data[2],
                  rot.data[4], rot.data[5], rot.data[6],
                  rot.data[8], rot.data[9], rot.data[10]);
        iitWorld = R * iitBody * R.transpose();
    }

    void calculateDerivedData() {
        orientation.normalise();
        transformMatrix.setOrientationAndPos(orientation, position);
        transformInertiaTensor(inverseInertiaTensorWorld, inverseInertiaTensor, transformMatrix);
    }

    void integrate(real dt) {
        if (!isAwakeFlag) return;
        lastFrameAcceleration = acceleration;
        lastFrameAcceleration.addScaledVector(forceAccum, inverseMass);
        Vector3 angularAcceleration = inverseInertiaTensorWorld.transform(torqueAccum);

        velocity.addScaledVector(lastFrameAcceleration, dt);
        rotation.addScaledVector(angularAcceleration, dt);
        velocity *= real_pow(linearDamping, dt);
        rotation *= real_pow(angularDamping, dt);
        position.addScaledVector(velocity, dt);
        orientation.addScaledVector(rotation, dt);

        calculateDerivedData();
        clearAccumulators();

        if (canSleepFlag) {
            real currentMotion = velocity.scalarProduct(velocity) + rotation.scalarProduct(rotation);
            real bias = real_pow((real)0.5, dt);
            motion = bias * motion + (1 - bias) * currentMotion;
            if (motion < sleepEpsilon) setAwake(false);
            else if (motion > 10 * sleepEpsilon) motion = 10 * sleepEpsilon;
        }
    }

    // -- force / torque accumulation --
    void addForce(const Vector3& f) { forceAccum += f; isAwakeFlag = true; }
    void addTorque(const Vector3& t) { torqueAccum += t; isAwakeFlag = true; }
    void addForceAtPoint(const Vector3& force, const Vector3& worldPoint) {
        Vector3 pt = worldPoint - position;
        forceAccum += force;
        torqueAccum += pt % force;
        isAwakeFlag = true;
    }
    void addForceAtBodyPoint(const Vector3& force, const Vector3& bodyPoint) {
        addForceAtPoint(force, getPointInWorldSpace(bodyPoint));
    }
    void clearAccumulators() { forceAccum.clear(); torqueAccum.clear(); }

    // -- transforms --
    Vector3 getPointInWorldSpace(const Vector3& p) const { return transformMatrix.transform(p); }
    Vector3 getPointInLocalSpace(const Vector3& p) const { return transformMatrix.transformInverse(p); }
    Vector3 getDirectionInWorldSpace(const Vector3& d) const { return transformMatrix.transformDirection(d); }
    Vector3 getDirectionInLocalSpace(const Vector3& d) const { return transformMatrix.transformInverseDirection(d); }

    // -- mass / inertia --
    void setMass(real m) { assert(m != 0); inverseMass = ((real)1) / m; }
    real getMass() const { return inverseMass == 0 ? REAL_MAX : ((real)1) / inverseMass; }
    void setInverseMass(real im) { inverseMass = im; }
    real getInverseMass() const { return inverseMass; }
    bool hasFiniteMass() const { return inverseMass > 0; }

    void setInertiaTensor(const Matrix3& inertiaTensor) { inverseInertiaTensor.setInverse(inertiaTensor); }
    void setInverseInertiaTensor(const Matrix3& iit) { inverseInertiaTensor = iit; }
    Matrix3 getInverseInertiaTensor() const { return inverseInertiaTensor; }
    Matrix3 getInverseInertiaTensorWorld() const { return inverseInertiaTensorWorld; }
    void getInverseInertiaTensorWorld(Matrix3* out) const { *out = inverseInertiaTensorWorld; }

    void setDamping(real lin, real ang) { linearDamping = lin; angularDamping = ang; }
    void setLinearDamping(real d) { linearDamping = d; }
    void setAngularDamping(real d) { angularDamping = d; }

    // -- state --
    void setPosition(const Vector3& p) { position = p; }
    void setPosition(real x, real y, real z) { position = Vector3(x, y, z); }
    Vector3 getPosition() const { return position; }

    void setOrientation(const Quaternion& q) { orientation = q; orientation.normalise(); }
    void setOrientation(real r, real i, real j, real k) { orientation = Quaternion(r, i, j, k); orientation.normalise(); }
    Quaternion getOrientation() const { return orientation; }
    void getOrientation(Matrix3* m) const { m->setOrientation(orientation); }
    Matrix4 getTransform() const { return transformMatrix; }
    void getTransform(Matrix4* m) const { *m = transformMatrix; }

    void setVelocity(const Vector3& v) { velocity = v; }
    void setVelocity(real x, real y, real z) { velocity = Vector3(x, y, z); }
    Vector3 getVelocity() const { return velocity; }
    void addVelocity(const Vector3& dv) { velocity += dv; }

    void setRotation(const Vector3& r) { rotation = r; }
    void setRotation(real x, real y, real z) { rotation = Vector3(x, y, z); }
    Vector3 getRotation() const { return rotation; }
    void addRotation(const Vector3& dr) { rotation += dr; }

    void setAcceleration(const Vector3& a) { acceleration = a; }
    void setAcceleration(real x, real y, real z) { acceleration = Vector3(x, y, z); }
    Vector3 getAcceleration() const { return acceleration; }
    Vector3 getLastFrameAcceleration() const { return lastFrameAcceleration; }

    bool getAwake() const { return isAwakeFlag; }
    void setAwake(bool awake = true) {
        if (awake) { isAwakeFlag = true; motion = sleepEpsilon * 2; }
        else { isAwakeFlag = false; velocity.clear(); rotation.clear(); }
    }
    bool getCanSleep() const { return canSleepFlag; }
    void setCanSleep(bool canSleep = true) { canSleepFlag = canSleep; if (!canSleep && !isAwakeFlag) setAwake(true); }

    real getKineticEnergy() const {
        if (inverseMass <= 0) return 0;
        real lin = ((real)0.5) * getMass() * velocity.squareMagnitude();
        Matrix3 I; I.setInverse(inverseInertiaTensor);        // body inertia
        Vector3 L = I * rotation;
        real ang = ((real)0.5) * rotation.scalarProduct(L);
        return lin + ang;
    }
};

inline real RigidBody::sleepEpsilon = (real)0.3;

} // namespace phys
