// Rigid-body force generators (Millington Part III ch.11): gravity, a spring
// between two body points, aerodynamic surfaces (fixed and control-surface), and
// buoyancy. A registry applies them each frame, adding both force and torque.
#pragma once
#include "body.h"
#include <vector>

namespace phys {

class ForceGenerator {
public:
    virtual ~ForceGenerator() {}
    virtual void updateForce(RigidBody* body, real dt) = 0;
};

class ForceRegistry {
protected:
    struct Registration { RigidBody* body; ForceGenerator* fg; };
    std::vector<Registration> registrations;
public:
    void add(RigidBody* body, ForceGenerator* fg) { registrations.push_back({body, fg}); }
    void remove(RigidBody* body, ForceGenerator* fg) {
        for (auto it = registrations.begin(); it != registrations.end(); ++it)
            if (it->body == body && it->fg == fg) { registrations.erase(it); return; }
    }
    void clear() { registrations.clear(); }
    void updateForces(real dt) { for (auto& r : registrations) r.fg->updateForce(r.body, dt); }
};

class Gravity : public ForceGenerator {
    Vector3 gravity;
public:
    Gravity(const Vector3& g) : gravity(g) {}
    void updateForce(RigidBody* body, real) override {
        if (!body->hasFiniteMass()) return;
        body->addForce(gravity * body->getMass());
    }
};

// Spring joining a point on this body to a point on another body.
class Spring : public ForceGenerator {
    Vector3 connectionPoint;       // in local (body) coords
    Vector3 otherConnectionPoint;  // in the other body's local coords
    RigidBody* other;
    real springConstant, restLength;
public:
    Spring(const Vector3& localConnectionPt, RigidBody* other, const Vector3& otherConnectionPt,
           real springConstant, real restLength)
        : connectionPoint(localConnectionPt), otherConnectionPoint(otherConnectionPt),
          other(other), springConstant(springConstant), restLength(restLength) {}
    void updateForce(RigidBody* body, real) override {
        Vector3 lws = body->getPointInWorldSpace(connectionPoint);
        Vector3 ows = other->getPointInWorldSpace(otherConnectionPoint);
        Vector3 force = lws - ows;
        real magnitude = force.magnitude();
        if (magnitude == 0) return;
        magnitude = real_abs(magnitude - restLength) * springConstant;
        force.normalise(); force *= -magnitude;
        body->addForceAtPoint(force, lws);
    }
};

// Aerodynamic surface with a fixed aero tensor, relative to a wind speed.
class Aero : public ForceGenerator {
protected:
    Matrix3 tensor;                // relates air velocity to force, in body coords
    Vector3 position;              // where the surface is, in body coords
    const Vector3* windspeed;
public:
    Aero(const Matrix3& tensor, const Vector3& position, const Vector3* windspeed)
        : tensor(tensor), position(position), windspeed(windspeed) {}
    void updateForce(RigidBody* body, real dt) override { updateForceFromTensor(body, dt, tensor); }
protected:
    void updateForceFromTensor(RigidBody* body, real, const Matrix3& useTensor) {
        Vector3 velocity = body->getVelocity() + *windspeed;
        Vector3 bodyVel = body->getTransform().transformInverseDirection(velocity);
        Vector3 bodyForce = useTensor.transform(bodyVel);
        Vector3 force = body->getTransform().transformDirection(bodyForce);
        body->addForceAtBodyPoint(force, position);
    }
};

// Aero surface whose tensor is interpolated by a control input in [-1,1].
class AeroControl : public Aero {
protected:
    Matrix3 maxTensor, minTensor;
    real controlSetting = 0;
    Matrix3 getTensor() {
        if (controlSetting <= -1) return minTensor;
        if (controlSetting >= 1) return maxTensor;
        if (controlSetting < 0) return Matrix3::linearInterpolate(minTensor, tensor, controlSetting + 1);
        if (controlSetting > 0) return Matrix3::linearInterpolate(tensor, maxTensor, controlSetting);
        return tensor;
    }
public:
    AeroControl(const Matrix3& base, const Matrix3& min, const Matrix3& max, const Vector3& position, const Vector3* windspeed)
        : Aero(base, position, windspeed), maxTensor(max), minTensor(min) {}
    void setControl(real value) { controlSetting = value; }
    void updateForce(RigidBody* body, real dt) override { updateForceFromTensor(body, dt, getTensor()); }
};

class Buoyancy : public ForceGenerator {
    Vector3 centreOfBuoyancy;      // body coords
    real maxDepth, volume, waterHeight, liquidDensity;
public:
    Buoyancy(const Vector3& cOfB, real maxDepth, real volume, real waterHeight, real liquidDensity = 1000)
        : centreOfBuoyancy(cOfB), maxDepth(maxDepth), volume(volume),
          waterHeight(waterHeight), liquidDensity(liquidDensity) {}
    void updateForce(RigidBody* body, real) override {
        Vector3 pointInWorld = body->getPointInWorldSpace(centreOfBuoyancy);
        real depth = pointInWorld.y;
        if (depth >= waterHeight + maxDepth) return;
        Vector3 force(0, 0, 0);
        if (depth <= waterHeight - maxDepth) force.y = liquidDensity * volume;
        else force.y = liquidDensity * volume * (depth - maxDepth - waterHeight) / (-2 * maxDepth);
        body->addForceAtPoint(force, pointInWorld);
    }
};

} // namespace phys
