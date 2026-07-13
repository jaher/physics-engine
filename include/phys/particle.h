// A Particle: the simplest simulated object — position, velocity and mass, with
// no orientation (Millington, Part I ch.3-4). Uses semi-implicit Euler with a
// drag/damping term for numerical stability.
#pragma once
#include "core.h"
#include <cassert>

namespace phys {

class Particle {
protected:
    real inverseMass = 0;          // 0 ⇒ infinite mass (immovable)
    real damping = (real)0.995;
    Vector3 position;
    Vector3 velocity;
    Vector3 acceleration;          // constant (e.g. gravity) applied every frame
    Vector3 forceAccum;            // accumulated force, cleared each integrate
public:
    void integrate(real dt) {
        if (inverseMass <= 0) return;
        assert(dt > 0);
        position.addScaledVector(velocity, dt);
        Vector3 resultingAcc = acceleration;
        resultingAcc.addScaledVector(forceAccum, inverseMass);
        velocity.addScaledVector(resultingAcc, dt);
        velocity *= real_pow(damping, dt);
        clearAccumulator();
    }
    void addForce(const Vector3& f) { forceAccum += f; }
    void clearAccumulator() { forceAccum.clear(); }
    Vector3 getForceAccum() const { return forceAccum; }

    void setMass(real m) { assert(m != 0); inverseMass = ((real)1) / m; }
    real getMass() const { return inverseMass == 0 ? REAL_MAX : ((real)1) / inverseMass; }
    void setInverseMass(real im) { inverseMass = im; }
    real getInverseMass() const { return inverseMass; }
    bool hasFiniteMass() const { return inverseMass > 0; }

    void setDamping(real d) { damping = d; }
    real getDamping() const { return damping; }

    void setPosition(const Vector3& p) { position = p; }
    void setPosition(real x, real y, real z) { position = Vector3(x, y, z); }
    Vector3 getPosition() const { return position; }

    void setVelocity(const Vector3& v) { velocity = v; }
    void setVelocity(real x, real y, real z) { velocity = Vector3(x, y, z); }
    Vector3 getVelocity() const { return velocity; }

    void setAcceleration(const Vector3& a) { acceleration = a; }
    void setAcceleration(real x, real y, real z) { acceleration = Vector3(x, y, z); }
    Vector3 getAcceleration() const { return acceleration; }

    real getKineticEnergy() const {
        if (inverseMass <= 0) return 0;
        return ((real)0.5) * getMass() * velocity.squareMagnitude();
    }
};

} // namespace phys
