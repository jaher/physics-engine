// Rigid-body contacts and the collision-resolution pipeline (Millington Part V
// ch.14-15). Each contact resolves in a local basis (normal = x): velocity is
// resolved by impulses (with optional anisotropic friction) and interpenetration
// by nonlinear projection that splits the correction between linear and angular
// motion by each body's inertia along the normal.
#pragma once
#include "body.h"

namespace phys {

class Contact {
    friend class ContactResolver;
public:
    RigidBody* body[2] = {nullptr, nullptr};
    real friction = 0;
    real restitution = 0;
    Vector3 contactPoint;
    Vector3 contactNormal;
    real penetration = 0;

    void setBodyData(RigidBody* one, RigidBody* two, real fric, real rest) {
        body[0] = one; body[1] = two; friction = fric; restitution = rest;
    }
protected:
    Matrix3 contactToWorld;                 // basis: columns are normal + two tangents
    Vector3 contactVelocity;                // closing velocity in contact coords
    real desiredDeltaVelocity = 0;
    Vector3 relativeContactPosition[2];

    void calculateInternals(real dt) {
        if (!body[0]) swapBodies();
        calculateContactBasis();
        relativeContactPosition[0] = contactPoint - body[0]->getPosition();
        if (body[1]) relativeContactPosition[1] = contactPoint - body[1]->getPosition();
        contactVelocity = calculateLocalVelocity(0, dt);
        if (body[1]) contactVelocity -= calculateLocalVelocity(1, dt);
        calculateDesiredDeltaVelocity(dt);
    }
    void swapBodies() { contactNormal *= -1; RigidBody* t = body[0]; body[0] = body[1]; body[1] = t; }

    void matchAwakeState() {
        if (!body[1]) return;
        bool a0 = body[0]->getAwake(), a1 = body[1]->getAwake();
        if (a0 ^ a1) { if (a0) body[1]->setAwake(true); else body[0]->setAwake(true); }
    }

    void calculateContactBasis() {
        Vector3 n = contactNormal, tangent[2];
        if (real_abs(n.x) > real_abs(n.y)) {
            const real s = ((real)1) / real_sqrt(n.z * n.z + n.x * n.x);
            tangent[0] = Vector3(n.z * s, 0, -n.x * s);
            tangent[1] = Vector3(n.y * tangent[0].x, n.z * tangent[0].x - n.x * tangent[0].z, -n.y * tangent[0].x);
        } else {
            const real s = ((real)1) / real_sqrt(n.z * n.z + n.y * n.y);
            tangent[0] = Vector3(0, -n.z * s, n.y * s);
            tangent[1] = Vector3(n.y * tangent[0].z - n.z * tangent[0].y, -n.x * tangent[0].z, n.x * tangent[0].y);
        }
        contactToWorld.setComponents(n, tangent[0], tangent[1]);
    }

    Vector3 calculateLocalVelocity(unsigned i, real dt) {
        RigidBody* b = body[i];
        Vector3 velocity = b->getRotation() % relativeContactPosition[i];
        velocity += b->getVelocity();
        Vector3 contactVel = contactToWorld.transformTranspose(velocity);
        Vector3 accVelocity = b->getLastFrameAcceleration() * dt;
        accVelocity = contactToWorld.transformTranspose(accVelocity);
        accVelocity.x = 0;                  // ignore the normal component of acceleration
        contactVel += accVelocity;
        return contactVel;
    }
    void calculateDesiredDeltaVelocity(real dt) {
        const static real velocityLimit = (real)0.25;
        real velocityFromAcc = 0;
        if (body[0]->getAwake()) velocityFromAcc += (body[0]->getLastFrameAcceleration() * dt) * contactNormal;
        if (body[1] && body[1]->getAwake()) velocityFromAcc -= (body[1]->getLastFrameAcceleration() * dt) * contactNormal;
        real thisRestitution = restitution;
        if (real_abs(contactVelocity.x) < velocityLimit) thisRestitution = 0;
        desiredDeltaVelocity = -contactVelocity.x - thisRestitution * (contactVelocity.x - velocityFromAcc);
    }

    Vector3 calculateFrictionlessImpulse(Matrix3* iit) {
        Vector3 deltaVelWorld = relativeContactPosition[0] % contactNormal;
        deltaVelWorld = iit[0].transform(deltaVelWorld);
        deltaVelWorld = deltaVelWorld % relativeContactPosition[0];
        real deltaVelocity = deltaVelWorld * contactNormal;
        deltaVelocity += body[0]->getInverseMass();
        if (body[1]) {
            Vector3 dvw = relativeContactPosition[1] % contactNormal;
            dvw = iit[1].transform(dvw);
            dvw = dvw % relativeContactPosition[1];
            deltaVelocity += dvw * contactNormal;
            deltaVelocity += body[1]->getInverseMass();
        }
        return Vector3(desiredDeltaVelocity / deltaVelocity, 0, 0);
    }
    Vector3 calculateFrictionImpulse(Matrix3* iit) {
        real inverseMass = body[0]->getInverseMass();
        Matrix3 impulseToTorque; impulseToTorque.setSkewSymmetric(relativeContactPosition[0]);
        Matrix3 deltaVelWorld = impulseToTorque;
        deltaVelWorld *= iit[0];
        deltaVelWorld *= impulseToTorque;
        deltaVelWorld *= -1;
        if (body[1]) {
            Matrix3 itt2; itt2.setSkewSymmetric(relativeContactPosition[1]);
            Matrix3 dvw2 = itt2; dvw2 *= iit[1]; dvw2 *= itt2; dvw2 *= -1;
            deltaVelWorld += dvw2;
            inverseMass += body[1]->getInverseMass();
        }
        Matrix3 deltaVelocity = contactToWorld.transpose();
        deltaVelocity *= deltaVelWorld;
        deltaVelocity *= contactToWorld;
        deltaVelocity.data[0] += inverseMass;
        deltaVelocity.data[4] += inverseMass;
        deltaVelocity.data[8] += inverseMass;
        Matrix3 impulseMatrix = deltaVelocity.inverse();
        Vector3 velKill(desiredDeltaVelocity, -contactVelocity.y, -contactVelocity.z);
        Vector3 impulseContact = impulseMatrix.transform(velKill);
        real planarImpulse = real_sqrt(impulseContact.y * impulseContact.y + impulseContact.z * impulseContact.z);
        if (planarImpulse > impulseContact.x * friction) {               // exceeded static friction → dynamic
            impulseContact.y /= planarImpulse;
            impulseContact.z /= planarImpulse;
            impulseContact.x = deltaVelocity.data[0]
                + deltaVelocity.data[1] * friction * impulseContact.y
                + deltaVelocity.data[2] * friction * impulseContact.z;
            impulseContact.x = desiredDeltaVelocity / impulseContact.x;
            impulseContact.y *= friction * impulseContact.x;
            impulseContact.z *= friction * impulseContact.x;
        }
        return impulseContact;
    }

    void applyVelocityChange(Vector3 velocityChange[2], Vector3 rotationChange[2]) {
        Matrix3 iit[2];
        body[0]->getInverseInertiaTensorWorld(&iit[0]);
        if (body[1]) body[1]->getInverseInertiaTensorWorld(&iit[1]);
        Vector3 impulseContact = (friction == 0) ? calculateFrictionlessImpulse(iit) : calculateFrictionImpulse(iit);
        Vector3 impulse = contactToWorld.transform(impulseContact);

        Vector3 impulsiveTorque = relativeContactPosition[0] % impulse;
        rotationChange[0] = iit[0].transform(impulsiveTorque);
        velocityChange[0].clear();
        velocityChange[0].addScaledVector(impulse, body[0]->getInverseMass());
        body[0]->addVelocity(velocityChange[0]);
        body[0]->addRotation(rotationChange[0]);
        if (body[1]) {
            Vector3 t = impulse % relativeContactPosition[1];
            rotationChange[1] = iit[1].transform(t);
            velocityChange[1].clear();
            velocityChange[1].addScaledVector(impulse, -body[1]->getInverseMass());
            body[1]->addVelocity(velocityChange[1]);
            body[1]->addRotation(rotationChange[1]);
        }
    }

    void applyPositionChange(Vector3 linearChange[2], Vector3 angularChange[2], real penetration) {
        const real angularLimit = (real)0.2;
        real angularMove[2], linearMove[2], angularInertia[2], linearInertia[2];
        real totalInertia = 0;
        for (unsigned i = 0; i < 2; i++) if (body[i]) {
            Matrix3 iit; body[i]->getInverseInertiaTensorWorld(&iit);
            Vector3 angularInertiaWorld = relativeContactPosition[i] % contactNormal;
            angularInertiaWorld = iit.transform(angularInertiaWorld);
            angularInertiaWorld = angularInertiaWorld % relativeContactPosition[i];
            angularInertia[i] = angularInertiaWorld * contactNormal;
            linearInertia[i] = body[i]->getInverseMass();
            totalInertia += linearInertia[i] + angularInertia[i];
        }
        for (unsigned i = 0; i < 2; i++) if (body[i]) {
            real sign = (i == 0) ? 1 : -1;
            angularMove[i] = sign * penetration * (angularInertia[i] / totalInertia);
            linearMove[i] = sign * penetration * (linearInertia[i] / totalInertia);

            Vector3 projection = relativeContactPosition[i];
            projection.addScaledVector(contactNormal, -relativeContactPosition[i].scalarProduct(contactNormal));
            real maxMagnitude = angularLimit * projection.magnitude();
            if (angularMove[i] < -maxMagnitude) { real total = angularMove[i] + linearMove[i]; angularMove[i] = -maxMagnitude; linearMove[i] = total - angularMove[i]; }
            else if (angularMove[i] > maxMagnitude) { real total = angularMove[i] + linearMove[i]; angularMove[i] = maxMagnitude; linearMove[i] = total - angularMove[i]; }

            if (angularMove[i] == 0) angularChange[i].clear();
            else {
                Vector3 targetAngularDirection = relativeContactPosition[i].vectorProduct(contactNormal);
                Matrix3 iit; body[i]->getInverseInertiaTensorWorld(&iit);
                angularChange[i] = iit.transform(targetAngularDirection) * (angularMove[i] / angularInertia[i]);
            }
            linearChange[i] = contactNormal * linearMove[i];

            Vector3 pos = body[i]->getPosition();
            pos.addScaledVector(contactNormal, linearMove[i]);
            body[i]->setPosition(pos);
            Quaternion q = body[i]->getOrientation();
            q.addScaledVector(angularChange[i], ((real)1));
            body[i]->setOrientation(q);
            body[i]->calculateDerivedData();
        }
    }
};

class ContactResolver {
protected:
    unsigned velocityIterations, positionIterations;
    real velocityEpsilon, positionEpsilon;
public:
    unsigned velocityIterationsUsed = 0, positionIterationsUsed = 0;

    ContactResolver(unsigned iterations, real velEps = (real)0.01, real posEps = (real)0.01) { setIterations(iterations); setEpsilon(velEps, posEps); }
    ContactResolver(unsigned vIt, unsigned pIt, real velEps = (real)0.01, real posEps = (real)0.01) { setIterations(vIt, pIt); setEpsilon(velEps, posEps); }
    void setIterations(unsigned it) { setIterations(it, it); }
    void setIterations(unsigned v, unsigned p) { velocityIterations = v; positionIterations = p; }
    void setEpsilon(real v, real p) { velocityEpsilon = v; positionEpsilon = p; }

    void resolveContacts(Contact* contacts, unsigned numContacts, real dt) {
        if (numContacts == 0) return;
        prepareContacts(contacts, numContacts, dt);
        adjustPositions(contacts, numContacts, dt);
        adjustVelocities(contacts, numContacts, dt);
    }
protected:
    void prepareContacts(Contact* contacts, unsigned numContacts, real dt) {
        for (unsigned i = 0; i < numContacts; i++) contacts[i].calculateInternals(dt);
    }
    void adjustVelocities(Contact* contacts, unsigned numContacts, real dt) {
        Vector3 velocityChange[2], rotationChange[2];
        velocityIterationsUsed = 0;
        while (velocityIterationsUsed < velocityIterations) {
            real max = velocityEpsilon; unsigned index = numContacts;
            for (unsigned i = 0; i < numContacts; i++)
                if (contacts[i].desiredDeltaVelocity > max) { max = contacts[i].desiredDeltaVelocity; index = i; }
            if (index == numContacts) break;
            contacts[index].matchAwakeState();
            contacts[index].applyVelocityChange(velocityChange, rotationChange);
            for (unsigned i = 0; i < numContacts; i++)
                for (unsigned b = 0; b < 2; b++) if (contacts[i].body[b])
                    for (unsigned d = 0; d < 2; d++)
                        if (contacts[i].body[b] == contacts[index].body[d]) {
                            Vector3 deltaVel = velocityChange[d] + rotationChange[d].vectorProduct(contacts[i].relativeContactPosition[b]);
                            contacts[i].contactVelocity += contacts[i].contactToWorld.transformTranspose(deltaVel) * (b ? -1 : 1);
                            contacts[i].calculateDesiredDeltaVelocity(dt);
                        }
            velocityIterationsUsed++;
        }
    }
    void adjustPositions(Contact* contacts, unsigned numContacts, real dt) {
        Vector3 linearChange[2], angularChange[2];
        positionIterationsUsed = 0;
        while (positionIterationsUsed < positionIterations) {
            real max = positionEpsilon; unsigned index = numContacts;
            for (unsigned i = 0; i < numContacts; i++)
                if (contacts[i].penetration > max) { max = contacts[i].penetration; index = i; }
            if (index == numContacts) break;
            contacts[index].matchAwakeState();
            contacts[index].applyPositionChange(linearChange, angularChange, max);
            for (unsigned i = 0; i < numContacts; i++)
                for (unsigned b = 0; b < 2; b++) if (contacts[i].body[b])
                    for (unsigned d = 0; d < 2; d++)
                        if (contacts[i].body[b] == contacts[index].body[d]) {
                            Vector3 deltaPosition = linearChange[d] + angularChange[d].vectorProduct(contacts[i].relativeContactPosition[b]);
                            contacts[i].penetration += deltaPosition.scalarProduct(contacts[i].contactNormal) * (b ? 1 : -1);
                        }
            positionIterationsUsed++;
        }
    }
};

// A source of contacts for the world to collect each frame.
class ContactGenerator {
public:
    virtual ~ContactGenerator() {}
    virtual unsigned addContact(Contact* contact, unsigned limit) const = 0;
};

} // namespace phys
