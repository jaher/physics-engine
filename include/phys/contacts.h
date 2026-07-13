// Rigid-body contacts and the collision-resolution pipeline (Millington Part V
// ch.14-15). Each contact resolves in a local basis (normal = x): velocity is
// resolved by impulses (with optional anisotropic friction) and interpenetration
// by nonlinear projection that splits the correction between linear and angular
// motion by each body's inertia along the normal.
#pragma once
#include "body.h"
#include "material.h"
#include <vector>

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
    // Set bodies + all four combined coefficients from a Material (see material.h).
    void setMaterial(RigidBody* one, RigidBody* two, const Material& m) {
        body[0] = one; body[1] = two;
        friction = m.friction; restitution = m.restitution;
        rollingFriction = m.rollingFriction; spinFriction = m.spinFriction;
    }
    // Enable rolling / spinning (torsional) friction on this contact.
    void setRollingFriction(real roll, real spin = 0) { rollingFriction = roll; spinFriction = spin; }
protected:
    Matrix3 contactToWorld;                 // basis: columns are normal + two tangents
    Vector3 contactVelocity;                // closing velocity in contact coords
    real desiredDeltaVelocity = 0;
    Vector3 relativeContactPosition[2];

    // -- persistent-manifold warm-starting & rolling/spinning friction --
    real rollingFriction = 0;               // 0 = disabled (torque opposing roll)
    real spinFriction = 0;                  // 0 = disabled (torque opposing spin about normal)
    real wsNormalImpulse = 0;               // previous frame's normal impulse (warm-start seed)
    Vector3 wsTangentImpulse;               // previous frame's friction impulse (world space)
    real accNormalImpulse = 0;              // normal impulse accumulated this frame
    Vector3 accTangentImpulse;              // friction impulse accumulated this frame (world)

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
        accNormalImpulse += impulseContact.x;
        accTangentImpulse += contactToWorld.transform(Vector3(0, impulseContact.y, impulseContact.z));
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

    // Effective inverse mass along the contact normal — the scalar the frictionless
    // solver divides by. Used to clamp warm-start impulses so they never inject energy.
    real normalDenominator() const {
        Matrix3 iit; body[0]->getInverseInertiaTensorWorld(&iit);
        Vector3 d = relativeContactPosition[0] % contactNormal;
        d = iit.transform(d); d = d % relativeContactPosition[0];
        real denom = d * contactNormal + body[0]->getInverseMass();
        if (body[1]) {
            Matrix3 iit2; body[1]->getInverseInertiaTensorWorld(&iit2);
            Vector3 d2 = relativeContactPosition[1] % contactNormal;
            d2 = iit2.transform(d2); d2 = d2 % relativeContactPosition[1];
            denom += d2 * contactNormal + body[1]->getInverseMass();
        }
        return denom;
    }

    // Apply an impulse given in contact coordinates (x = normal, y,z = tangent),
    // recording it in this frame's accumulated impulses and returning the per-body
    // velocity / rotation changes so the caller can propagate them to other contacts.
    void applyContactImpulse(const Vector3& impulseContact, Vector3 velocityChange[2], Vector3 rotationChange[2]) {
        Matrix3 iit[2];
        body[0]->getInverseInertiaTensorWorld(&iit[0]);
        if (body[1]) body[1]->getInverseInertiaTensorWorld(&iit[1]);
        Vector3 impulse = contactToWorld.transform(impulseContact);
        Vector3 impulsiveTorque = relativeContactPosition[0] % impulse;
        rotationChange[0] = iit[0].transform(impulsiveTorque);
        velocityChange[0].clear();
        velocityChange[0].addScaledVector(impulse, body[0]->getInverseMass());
        body[0]->addVelocity(velocityChange[0]);
        body[0]->addRotation(rotationChange[0]);
        velocityChange[1].clear(); rotationChange[1].clear();
        if (body[1]) {
            Vector3 t = impulse % relativeContactPosition[1];
            rotationChange[1] = iit[1].transform(t);
            velocityChange[1].addScaledVector(impulse, -body[1]->getInverseMass());
            body[1]->addVelocity(velocityChange[1]);
            body[1]->addRotation(rotationChange[1]);
        }
        accNormalImpulse += impulseContact.x;
        accTangentImpulse += contactToWorld.transform(Vector3(0, impulseContact.y, impulseContact.z));
    }

    // Rolling and spinning (torsional) friction. Dissipates the tangential (rolling)
    // and normal (spinning) components of the relative angular velocity, each
    // Coulomb-limited by the accumulated normal impulse. Applied purely as an
    // angular impulse, so it lets a rolling sphere slow to rest without a linear kick.
    void applyRollingFriction(real /*dt*/) {
        if ((rollingFriction <= 0 && spinFriction <= 0) || accNormalImpulse <= 0) return;
        Matrix3 iit0; body[0]->getInverseInertiaTensorWorld(&iit0);
        Matrix3 iit1; bool has1 = body[1] != nullptr;
        if (has1) body[1]->getInverseInertiaTensorWorld(&iit1);
        Vector3 relOmega = body[0]->getRotation();
        if (has1) relOmega -= body[1]->getRotation();
        Vector3 n = contactNormal;
        real spin = relOmega * n;                    // component along the normal
        Vector3 rollOmega = relOmega - n * spin;     // tangential (rolling) component

        if (rollingFriction > 0) {
            real rmag = rollOmega.magnitude();
            if (rmag > (real)1e-9) {
                Vector3 axis = rollOmega * (((real)1) / rmag);
                real effI = axis * iit0.transform(axis);
                if (has1) effI += axis * iit1.transform(axis);
                if (effI > 0) {
                    real full = rmag / effI;                        // impulse to fully stop
                    real maxImp = rollingFriction * accNormalImpulse;
                    real lambda = full < maxImp ? full : maxImp;
                    Vector3 L = axis * (-lambda);
                    body[0]->addRotation(iit0.transform(L));
                    if (has1) body[1]->addRotation(iit1.transform(-L));
                }
            }
        }
        if (spinFriction > 0) {
            real smag = real_abs(spin);
            if (smag > (real)1e-9) {
                Vector3 axis = n * (spin > 0 ? (real)1 : (real)-1);
                real effI = axis * iit0.transform(axis);
                if (has1) effI += axis * iit1.transform(axis);
                if (effI > 0) {
                    real full = smag / effI;
                    real maxImp = spinFriction * accNormalImpulse;
                    real lambda = full < maxImp ? full : maxImp;
                    Vector3 L = axis * (-lambda);
                    body[0]->addRotation(iit0.transform(L));
                    if (has1) body[1]->addRotation(iit1.transform(-L));
                }
            }
        }
    }
};

// One persistent contact point in a manifold: its position in each body's local
// frame (for cross-frame matching) and the impulses accumulated last frame.
struct ManifoldPoint {
    Vector3 localA, localB;
    real normalImpulse = 0;
    Vector3 tangentImpulse;                 // world space
    bool seen = false;
};

// Up to four persistent contact points shared by one body-pair, kept across frames
// so the solver can warm-start from the previous frame's impulses (Bullet/PhysX style).
struct ContactManifold {
    RigidBody* body[2] = {nullptr, nullptr};
    ManifoldPoint pts[4];
    unsigned count = 0;
};

// Maps a contact (this frame) to its slot (manifold index, point index) in the cache.
struct ManifoldSlot { int mf = -1, pt = -1; };

class ContactResolver {
protected:
    unsigned velocityIterations, positionIterations;
    real velocityEpsilon, positionEpsilon;
    std::vector<ContactManifold> manifolds;     // persistent contact cache (warm-start)
    std::vector<ManifoldSlot> contactSlot;       // this-frame contact → cache slot map
    bool warmStarting = true;
public:
    unsigned velocityIterationsUsed = 0, positionIterationsUsed = 0;

    ContactResolver(unsigned iterations, real velEps = (real)0.01, real posEps = (real)0.01) { setIterations(iterations); setEpsilon(velEps, posEps); }
    ContactResolver(unsigned vIt, unsigned pIt, real velEps = (real)0.01, real posEps = (real)0.01) { setIterations(vIt, pIt); setEpsilon(velEps, posEps); }
    void setIterations(unsigned it) { setIterations(it, it); }
    void setIterations(unsigned v, unsigned p) { velocityIterations = v; positionIterations = p; }
    void setEpsilon(real v, real p) { velocityEpsilon = v; positionEpsilon = p; }

    // Warm-starting from the persistent manifold cache (on by default). Toggle off
    // for a byte-for-byte cold solve (used to benchmark the two against each other).
    bool getWarmStarting() const { return warmStarting; }
    void setWarmStarting(bool w) { warmStarting = w; }
    void clearManifolds() { manifolds.clear(); }

    void resolveContacts(Contact* contacts, unsigned numContacts, real dt) {
        if (numContacts == 0) return;
        prepareContacts(contacts, numContacts, dt);
        if (warmStarting) matchManifolds(contacts, numContacts);
        adjustPositions(contacts, numContacts, dt);
        if (warmStarting) applyWarmStart(contacts, numContacts, dt);
        adjustVelocities(contacts, numContacts, dt);
        applyRollingFriction(contacts, numContacts, dt);
        if (warmStarting) storeManifolds(contacts, numContacts);
    }
protected:
    void prepareContacts(Contact* contacts, unsigned numContacts, real dt) {
        for (unsigned i = 0; i < numContacts; i++) {
            contacts[i].calculateInternals(dt);
            contacts[i].wsNormalImpulse = 0;  contacts[i].wsTangentImpulse = Vector3();
            contacts[i].accNormalImpulse = 0; contacts[i].accTangentImpulse = Vector3();
        }
    }

    // Match this frame's contacts to persistent manifold points by body-pair and
    // body-local proximity, seeding each contact's warm-start impulse from the cache.
    void matchManifolds(Contact* contacts, unsigned numContacts) {
        const real matchTol2 = (real)0.0004;    // (2 cm)^2 in body-local space
        for (ContactManifold& m : manifolds)
            for (unsigned k = 0; k < m.count; k++) m.pts[k].seen = false;
        contactSlot.assign(numContacts, ManifoldSlot());
        for (unsigned i = 0; i < numContacts; i++) {
            Contact& c = contacts[i];
            if (!c.body[0]) continue;
            Vector3 la = c.body[0]->getPointInLocalSpace(c.contactPoint);
            Vector3 lb = c.body[1] ? c.body[1]->getPointInLocalSpace(c.contactPoint) : Vector3();
            int mi = -1;
            for (unsigned k = 0; k < manifolds.size(); k++)
                if (manifolds[k].body[0] == c.body[0] && manifolds[k].body[1] == c.body[1]) { mi = (int)k; break; }
            if (mi < 0) {
                ContactManifold m; m.body[0] = c.body[0]; m.body[1] = c.body[1];
                manifolds.push_back(m); mi = (int)manifolds.size() - 1;
            }
            ContactManifold& m = manifolds[mi];
            int pj = -1; real best = matchTol2;
            for (unsigned k = 0; k < m.count; k++) {
                if (m.pts[k].seen) continue;
                real d = (m.pts[k].localA - la).squareMagnitude();
                if (c.body[1]) d += (m.pts[k].localB - lb).squareMagnitude();
                if (d < best) { best = d; pj = (int)k; }
            }
            if (pj < 0) {
                if (m.count < 4) { pj = (int)m.count; m.pts[pj] = ManifoldPoint(); m.count++; }
                else pj = 0;                          // full & no match → reuse slot 0
            }
            ManifoldPoint& p = m.pts[pj];
            c.wsNormalImpulse = p.normalImpulse;      // 0 for a freshly created point
            c.wsTangentImpulse = p.tangentImpulse;
            p.localA = la; p.localB = lb; p.seen = true;
            contactSlot[i].mf = mi; contactSlot[i].pt = pj;
        }
    }

    // Seed the velocity solver with each matched contact's previous-frame impulses,
    // clamped so they can only slow an approaching / sliding contact (never inject
    // energy) and skipped for fast bounces. Propagates like a normal solver iteration.
    void applyWarmStart(Contact* contacts, unsigned numContacts, real dt) {
        const real velocityLimit = (real)0.25;
        Vector3 velocityChange[2], rotationChange[2];
        for (unsigned index = 0; index < numContacts; index++) {
            Contact& c = contacts[index];
            if (c.wsNormalImpulse <= 0 && c.wsTangentImpulse.squareMagnitude() <= 0) continue;
            real vN = c.contactVelocity.x;                       // <0 approaching, >0 separating
            if (vN >= 0 || vN <= -velocityLimit) continue;       // separating, or a fast bounce
            real denom = c.normalDenominator();
            if (denom <= 0) continue;
            real target = -vN / denom;                           // impulse that brings closing vel to 0
            real Pn = c.wsNormalImpulse; if (Pn > target) Pn = target; if (Pn < 0) Pn = 0;
            real py = 0, pz = 0;                                  // friction warm-start (contact tangent)
            real vy = c.contactVelocity.y, vz = c.contactVelocity.z;
            real vT = real_sqrt(vy * vy + vz * vz);
            if (Pn > 0 && vT > (real)1e-9 && c.friction > 0) {
                Vector3 wt = c.contactToWorld.transformTranspose(c.wsTangentImpulse);
                real prevMag = real_sqrt(wt.y * wt.y + wt.z * wt.z);
                real mag = vT / denom;                           // impulse that stops tangential slide
                real cone = c.friction * Pn;
                if (cone < mag) mag = cone;                      // Coulomb cone
                if (prevMag < mag) mag = prevMag;                // history-limited (warm-start seed)
                py = -vy / vT * mag; pz = -vz / vT * mag;
            }
            if (Pn <= 0 && py == 0 && pz == 0) continue;
            c.applyContactImpulse(Vector3(Pn, py, pz), velocityChange, rotationChange);
            for (unsigned i = 0; i < numContacts; i++)
                for (unsigned b = 0; b < 2; b++) if (contacts[i].body[b])
                    for (unsigned d = 0; d < 2; d++)
                        if (contacts[i].body[b] == contacts[index].body[d]) {
                            Vector3 deltaVel = velocityChange[d] + rotationChange[d].vectorProduct(contacts[i].relativeContactPosition[b]);
                            contacts[i].contactVelocity += contacts[i].contactToWorld.transformTranspose(deltaVel) * (b ? -1 : 1);
                            contacts[i].calculateDesiredDeltaVelocity(dt);
                        }
        }
    }

    // Write the frame's accumulated impulses back into the cache, then drop points /
    // manifolds that saw no contact this frame so the cache stays bounded and fresh.
    void storeManifolds(Contact* contacts, unsigned numContacts) {
        for (unsigned i = 0; i < numContacts; i++) {
            ManifoldSlot s = contactSlot[i];
            if (s.mf < 0) continue;
            ManifoldPoint& p = manifolds[s.mf].pts[s.pt];
            p.normalImpulse = contacts[i].accNormalImpulse;
            p.tangentImpulse = contacts[i].accTangentImpulse;
        }
        for (unsigned k = 0; k < manifolds.size(); ) {
            ContactManifold& m = manifolds[k];
            unsigned w = 0;
            for (unsigned j = 0; j < m.count; j++) if (m.pts[j].seen) m.pts[w++] = m.pts[j];
            m.count = w;
            if (w == 0) { manifolds[k] = manifolds.back(); manifolds.pop_back(); }
            else k++;
        }
    }

    void applyRollingFriction(Contact* contacts, unsigned numContacts, real dt) {
        for (unsigned i = 0; i < numContacts; i++) contacts[i].applyRollingFriction(dt);
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
