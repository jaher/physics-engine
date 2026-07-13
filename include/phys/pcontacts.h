// Particle contacts and their resolver (Millington Part II ch.7). A contact
// holds two particles (the second may be null for a scenery contact), a
// restitution, a normal and a penetration; resolution handles both the velocity
// (impulse) and the interpenetration (projection), with the resting-contact
// correction for accelerated approach.
#pragma once
#include "particle.h"

namespace phys {

class ParticleContact {
    friend class ParticleContactResolver;
public:
    Particle* particle[2] = {nullptr, nullptr};
    real restitution = 0;
    Vector3 contactNormal;
    real penetration = 0;
    Vector3 particleMovement[2];

    void resolve(real dt) { resolveVelocity(dt); resolveInterpenetration(dt); }

    real calculateSeparatingVelocity() const {
        Vector3 relativeVelocity = particle[0]->getVelocity();
        if (particle[1]) relativeVelocity -= particle[1]->getVelocity();
        return relativeVelocity * contactNormal;
    }
private:
    void resolveVelocity(real dt) {
        real separatingVelocity = calculateSeparatingVelocity();
        if (separatingVelocity > 0) return;               // separating or resting → no impulse
        real newSepVelocity = -separatingVelocity * restitution;

        // remove velocity built up by acceleration this frame (resting contacts)
        Vector3 accCausedVelocity = particle[0]->getAcceleration();
        if (particle[1]) accCausedVelocity -= particle[1]->getAcceleration();
        real accCausedSepVelocity = (accCausedVelocity * contactNormal) * dt;
        if (accCausedSepVelocity < 0) {
            newSepVelocity += restitution * accCausedSepVelocity;
            if (newSepVelocity < 0) newSepVelocity = 0;
        }
        real deltaVelocity = newSepVelocity - separatingVelocity;

        real totalInverseMass = particle[0]->getInverseMass();
        if (particle[1]) totalInverseMass += particle[1]->getInverseMass();
        if (totalInverseMass <= 0) return;

        real impulse = deltaVelocity / totalInverseMass;
        Vector3 impulsePerIMass = contactNormal * impulse;
        particle[0]->setVelocity(particle[0]->getVelocity() + impulsePerIMass * particle[0]->getInverseMass());
        if (particle[1])
            particle[1]->setVelocity(particle[1]->getVelocity() + impulsePerIMass * -particle[1]->getInverseMass());
    }
    void resolveInterpenetration(real) {
        if (penetration <= 0) { particleMovement[0].clear(); particleMovement[1].clear(); return; }
        real totalInverseMass = particle[0]->getInverseMass();
        if (particle[1]) totalInverseMass += particle[1]->getInverseMass();
        if (totalInverseMass <= 0) { particleMovement[0].clear(); particleMovement[1].clear(); return; }

        Vector3 movePerIMass = contactNormal * (penetration / totalInverseMass);
        particleMovement[0] = movePerIMass * particle[0]->getInverseMass();
        particleMovement[1] = particle[1] ? movePerIMass * -particle[1]->getInverseMass() : Vector3();
        particle[0]->setPosition(particle[0]->getPosition() + particleMovement[0]);
        if (particle[1]) particle[1]->setPosition(particle[1]->getPosition() + particleMovement[1]);
    }
};

class ParticleContactResolver {
protected:
    unsigned iterations;
    unsigned iterationsUsed = 0;
public:
    ParticleContactResolver(unsigned iterations) : iterations(iterations) {}
    void setIterations(unsigned it) { iterations = it; }

    void resolveContacts(ParticleContact* contactArray, unsigned numContacts, real dt) {
        iterationsUsed = 0;
        while (iterationsUsed < iterations) {
            real max = REAL_MAX; unsigned maxIndex = numContacts;
            for (unsigned i = 0; i < numContacts; i++) {
                real sepVel = contactArray[i].calculateSeparatingVelocity();
                if (sepVel < max && (sepVel < 0 || contactArray[i].penetration > 0)) { max = sepVel; maxIndex = i; }
            }
            if (maxIndex == numContacts) break;

            contactArray[maxIndex].resolve(dt);
            Vector3* move = contactArray[maxIndex].particleMovement;
            for (unsigned i = 0; i < numContacts; i++) {
                if (contactArray[i].particle[0] == contactArray[maxIndex].particle[0])
                    contactArray[i].penetration -= move[0] * contactArray[i].contactNormal;
                else if (contactArray[i].particle[0] == contactArray[maxIndex].particle[1])
                    contactArray[i].penetration -= move[1] * contactArray[i].contactNormal;
                if (contactArray[i].particle[1]) {
                    if (contactArray[i].particle[1] == contactArray[maxIndex].particle[0])
                        contactArray[i].penetration += move[0] * contactArray[i].contactNormal;
                    else if (contactArray[i].particle[1] == contactArray[maxIndex].particle[1])
                        contactArray[i].penetration += move[1] * contactArray[i].contactNormal;
                }
            }
            iterationsUsed++;
        }
    }
};

class ParticleContactGenerator {
public:
    virtual ~ParticleContactGenerator() {}
    virtual unsigned addContact(ParticleContact* contact, unsigned limit) const = 0;
};

} // namespace phys
