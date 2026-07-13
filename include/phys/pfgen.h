// Particle force generators (Millington Part II ch.5-6): gravity, drag, the
// spring family (spring, anchored spring, bungee, buoyancy, and a stiff "fake"
// spring solved analytically) plus a registry that applies them each frame.
#pragma once
#include "particle.h"
#include <vector>

namespace phys {

class ParticleForceGenerator {
public:
    virtual ~ParticleForceGenerator() {}
    virtual void updateForce(Particle* particle, real dt) = 0;
};

class ParticleForceRegistry {
protected:
    struct Registration { Particle* particle; ParticleForceGenerator* fg; };
    std::vector<Registration> registrations;
public:
    void add(Particle* particle, ParticleForceGenerator* fg) { registrations.push_back({particle, fg}); }
    void remove(Particle* particle, ParticleForceGenerator* fg) {
        for (auto it = registrations.begin(); it != registrations.end(); ++it)
            if (it->particle == particle && it->fg == fg) { registrations.erase(it); return; }
    }
    void clear() { registrations.clear(); }
    void updateForces(real dt) { for (auto& r : registrations) r.fg->updateForce(r.particle, dt); }
};

class ParticleGravity : public ParticleForceGenerator {
    Vector3 gravity;
public:
    ParticleGravity(const Vector3& g) : gravity(g) {}
    void updateForce(Particle* p, real) override {
        if (!p->hasFiniteMass()) return;
        p->addForce(gravity * p->getMass());
    }
};

class ParticleDrag : public ParticleForceGenerator {
    real k1, k2;                    // linear and quadratic drag coefficients
public:
    ParticleDrag(real k1, real k2) : k1(k1), k2(k2) {}
    void updateForce(Particle* p, real) override {
        Vector3 force = p->getVelocity();
        real speed = force.magnitude();
        if (speed == 0) return;
        real drag = k1 * speed + k2 * speed * speed;
        force.normalise();
        force *= -drag;
        p->addForce(force);
    }
};

class ParticleSpring : public ParticleForceGenerator {
    Particle* other; real springConstant, restLength;
public:
    ParticleSpring(Particle* other, real k, real rest) : other(other), springConstant(k), restLength(rest) {}
    void updateForce(Particle* p, real) override {
        Vector3 force = p->getPosition() - other->getPosition();
        real length = force.magnitude(); if (length == 0) return;
        real mag = springConstant * (length - restLength);
        force.normalise(); force *= -mag;
        p->addForce(force);
    }
};

class ParticleAnchoredSpring : public ParticleForceGenerator {
protected:
    Vector3* anchor; real springConstant, restLength;
public:
    ParticleAnchoredSpring() {}
    ParticleAnchoredSpring(Vector3* a, real k, real rest) : anchor(a), springConstant(k), restLength(rest) {}
    const Vector3* getAnchor() const { return anchor; }
    void init(Vector3* a, real k, real rest) { anchor = a; springConstant = k; restLength = rest; }
    void updateForce(Particle* p, real) override {
        Vector3 force = p->getPosition() - *anchor;
        real length = force.magnitude(); if (length == 0) return;
        real mag = (restLength - length) * springConstant;
        force.normalise(); force *= mag;
        p->addForce(force);
    }
};

class ParticleBungee : public ParticleForceGenerator {
    Particle* other; real springConstant, restLength;
public:
    ParticleBungee(Particle* other, real k, real rest) : other(other), springConstant(k), restLength(rest) {}
    void updateForce(Particle* p, real) override {
        Vector3 force = p->getPosition() - other->getPosition();
        real length = force.magnitude();
        if (length <= restLength) return;                 // bungees only pull
        real mag = springConstant * (length - restLength);
        force.normalise(); force *= -mag;
        p->addForce(force);
    }
};

class ParticleAnchoredBungee : public ParticleAnchoredSpring {
public:
    ParticleAnchoredBungee(Vector3* a, real k, real rest) : ParticleAnchoredSpring(a, k, rest) {}
    void updateForce(Particle* p, real) override {
        Vector3 force = p->getPosition() - *anchor;
        real length = force.magnitude();
        if (length <= restLength) return;
        real mag = springConstant * (length - restLength);
        force.normalise(); force *= -mag;
        p->addForce(force);
    }
};

class ParticleBuoyancy : public ParticleForceGenerator {
    real maxDepth, volume, waterHeight, liquidDensity;
public:
    ParticleBuoyancy(real maxDepth, real volume, real waterHeight, real liquidDensity = 1000)
        : maxDepth(maxDepth), volume(volume), waterHeight(waterHeight), liquidDensity(liquidDensity) {}
    void updateForce(Particle* p, real) override {
        real depth = p->getPosition().y;
        if (depth >= waterHeight + maxDepth) return;      // out of the water
        Vector3 force(0, 0, 0);
        if (depth <= waterHeight - maxDepth) { force.y = liquidDensity * volume; p->addForce(force); return; }
        force.y = liquidDensity * volume * (depth - maxDepth - waterHeight) / (-2 * maxDepth);
        p->addForce(force);
    }
};

// A very stiff spring solved analytically, avoiding the explosion an explicit
// integrator would suffer (assumes the anchor's inverse mass is negligible).
class ParticleFakeSpring : public ParticleForceGenerator {
    Vector3* anchor; real springConstant, damping;
public:
    ParticleFakeSpring(Vector3* a, real k, real d) : anchor(a), springConstant(k), damping(d) {}
    void updateForce(Particle* p, real dt) override {
        if (!p->hasFiniteMass()) return;
        Vector3 position = p->getPosition() - *anchor;
        real gamma = ((real)0.5) * real_sqrt(4 * springConstant - damping * damping);
        if (gamma == 0) return;
        Vector3 c = position * (damping / (2 * gamma)) + p->getVelocity() * (((real)1) / gamma);
        Vector3 target = position * real_cos(gamma * dt) + c * real_sin(gamma * dt);
        target *= real_exp(((real)-0.5) * dt * damping);
        Vector3 accel = (target - position) * (((real)1) / (dt * dt)) - p->getVelocity() * (((real)1) / dt);
        p->addForce(accel * p->getMass());
    }
};

} // namespace phys
