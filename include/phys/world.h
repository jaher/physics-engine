// The rigid-body world (Millington ch.15): owns bodies, their force generators
// and contact generators, and runs a whole physics frame — accumulate forces,
// integrate, generate contacts, resolve.
#pragma once
#include "fgen.h"
#include "contacts.h"
#include <vector>

namespace phys {

class World {
public:
    typedef std::vector<RigidBody*> RigidBodies;
    typedef std::vector<ContactGenerator*> ContactGenerators;
protected:
    RigidBodies bodies;
    ContactGenerators contactGenerators;
    ForceRegistry registry;
    ContactResolver resolver;
    Contact* contacts;
    unsigned maxContacts;
    bool calculateIterations;
public:
    World(unsigned maxContacts, unsigned iterations = 0)
        : resolver(iterations), maxContacts(maxContacts), calculateIterations(iterations == 0) {
        contacts = new Contact[maxContacts];
    }
    ~World() { delete[] contacts; }

    RigidBodies& getRigidBodies() { return bodies; }
    ContactGenerators& getContactGenerators() { return contactGenerators; }
    ForceRegistry& getForceRegistry() { return registry; }

    void startFrame() {
        for (RigidBody* b : bodies) { b->clearAccumulators(); b->calculateDerivedData(); }
    }
    unsigned generateContacts() {
        unsigned limit = maxContacts;
        Contact* nextContact = contacts;
        for (ContactGenerator* g : contactGenerators) {
            unsigned used = g->addContact(nextContact, limit);
            limit -= used; nextContact += used;
            if (limit <= 0) break;
        }
        return maxContacts - limit;
    }
    void integrate(real dt) { for (RigidBody* b : bodies) b->integrate(dt); }

    void runPhysics(real dt) {
        registry.updateForces(dt);
        integrate(dt);
        unsigned usedContacts = generateContacts();
        if (usedContacts) {
            if (calculateIterations) resolver.setIterations(usedContacts * 4, usedContacts * 4);
            resolver.resolveContacts(contacts, usedContacts, dt);
        }
    }
};

} // namespace phys
