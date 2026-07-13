// Hard constraints between particles expressed as contact generators
// (Millington Part II ch.7): cables (max length, bouncy), rods (fixed length),
// and their anchored-to-a-point variants.
#pragma once
#include "pcontacts.h"

namespace phys {

class ParticleLink : public ParticleContactGenerator {
public:
    Particle* particle[2] = {nullptr, nullptr};
protected:
    real currentLength() const { return (particle[0]->getPosition() - particle[1]->getPosition()).magnitude(); }
};

class ParticleCable : public ParticleLink {
public:
    real maxLength = 0, restitution = 0;
    unsigned addContact(ParticleContact* contact, unsigned limit) const override {
        real length = currentLength();
        if (length < maxLength) return 0;
        contact->particle[0] = particle[0]; contact->particle[1] = particle[1];
        Vector3 normal = particle[1]->getPosition() - particle[0]->getPosition();
        normal.normalise();
        contact->contactNormal = normal;
        contact->penetration = length - maxLength;
        contact->restitution = restitution;
        return 1;
    }
};

class ParticleRod : public ParticleLink {
public:
    real length = 0;
    unsigned addContact(ParticleContact* contact, unsigned limit) const override {
        real currentLen = currentLength();
        if (currentLen == length) return 0;
        contact->particle[0] = particle[0]; contact->particle[1] = particle[1];
        Vector3 normal = particle[1]->getPosition() - particle[0]->getPosition();
        normal.normalise();
        if (currentLen > length) { contact->contactNormal = normal; contact->penetration = currentLen - length; }
        else { contact->contactNormal = normal * -1; contact->penetration = length - currentLen; }
        contact->restitution = 0;                         // rods never bounce
        return 1;
    }
};

class ParticleConstraint : public ParticleContactGenerator {
public:
    Particle* particle = nullptr;
    Vector3 anchor;
protected:
    real currentLength() const { return (particle->getPosition() - anchor).magnitude(); }
};

class ParticleCableConstraint : public ParticleConstraint {
public:
    real maxLength = 0, restitution = 0;
    unsigned addContact(ParticleContact* contact, unsigned limit) const override {
        real length = currentLength();
        if (length < maxLength) return 0;
        contact->particle[0] = particle; contact->particle[1] = nullptr;
        Vector3 normal = anchor - particle->getPosition();
        normal.normalise();
        contact->contactNormal = normal;
        contact->penetration = length - maxLength;
        contact->restitution = restitution;
        return 1;
    }
};

class ParticleRodConstraint : public ParticleConstraint {
public:
    real length = 0;
    unsigned addContact(ParticleContact* contact, unsigned limit) const override {
        real currentLen = currentLength();
        if (currentLen == length) return 0;
        contact->particle[0] = particle; contact->particle[1] = nullptr;
        Vector3 normal = anchor - particle->getPosition();
        normal.normalise();
        if (currentLen > length) { contact->contactNormal = normal; contact->penetration = currentLen - length; }
        else { contact->contactNormal = normal * -1; contact->penetration = length - currentLen; }
        contact->restitution = 0;
        return 1;
    }
};

} // namespace phys
