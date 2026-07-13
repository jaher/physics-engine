// The mass-aggregate particle world (Millington Part II ch.8): keeps particles,
// their force generators and contact generators, and runs a full physics frame
// (forces → integrate → generate contacts → resolve).
#pragma once
#include "pfgen.h"
#include "plinks.h"
#include <vector>

namespace phys {

class ParticleWorld {
public:
    typedef std::vector<Particle*> Particles;
    typedef std::vector<ParticleContactGenerator*> ContactGenerators;
protected:
    Particles particles;
    ContactGenerators contactGenerators;
    ParticleForceRegistry registry;
    ParticleContactResolver resolver;
    ParticleContact* contacts;
    unsigned maxContacts;
    bool calculateIterations;
public:
    ParticleWorld(unsigned maxContacts, unsigned iterations = 0)
        : resolver(iterations), maxContacts(maxContacts),
          calculateIterations(iterations == 0) {
        contacts = new ParticleContact[maxContacts];
    }
    ~ParticleWorld() { delete[] contacts; }

    Particles& getParticles() { return particles; }
    ContactGenerators& getContactGenerators() { return contactGenerators; }
    ParticleForceRegistry& getForceRegistry() { return registry; }

    void startFrame() { for (Particle* p : particles) p->clearAccumulator(); }

    unsigned generateContacts() {
        unsigned limit = maxContacts;
        ParticleContact* nextContact = contacts;
        for (ParticleContactGenerator* g : contactGenerators) {
            unsigned used = g->addContact(nextContact, limit);
            limit -= used; nextContact += used;
            if (limit <= 0) break;
        }
        return maxContacts - limit;
    }
    void integrate(real dt) { for (Particle* p : particles) p->integrate(dt); }

    void runPhysics(real dt) {
        registry.updateForces(dt);
        integrate(dt);
        unsigned usedContacts = generateContacts();
        if (usedContacts) {
            if (calculateIterations) resolver.setIterations(usedContacts * 2);
            resolver.resolveContacts(contacts, usedContacts, dt);
        }
    }
};

// Scenery contact generator: keeps particles above the y = 0 ground plane.
class GroundContacts : public ParticleContactGenerator {
    ParticleWorld::Particles* particles = nullptr;
public:
    real restitution = (real)0.2;
    void init(ParticleWorld::Particles* p) { particles = p; }
    unsigned addContact(ParticleContact* contact, unsigned limit) const override {
        unsigned count = 0;
        for (Particle* p : *particles) {
            real y = p->getPosition().y;
            if (y < 0) {
                contact->contactNormal = Vector3::UP;
                contact->particle[0] = p; contact->particle[1] = nullptr;
                contact->penetration = -y;
                contact->restitution = restitution;
                contact++; count++;
                if (count >= limit) return count;
            }
        }
        return count;
    }
};

} // namespace phys
