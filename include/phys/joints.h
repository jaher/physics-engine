// A joint constraining a point on one body to coincide with a point on another
// (Millington ch.15). Implemented as a contact generator so the standard contact
// resolver keeps the two anchor points together.
#pragma once
#include "contacts.h"

namespace phys {

class Joint : public ContactGenerator {
public:
    RigidBody* body[2] = {nullptr, nullptr};
    Vector3 position[2];               // anchor in each body's local coords
    real error = 0;                    // permitted slack

    void set(RigidBody* a, const Vector3& aPos, RigidBody* b, const Vector3& bPos, real err) {
        body[0] = a; body[1] = b; position[0] = aPos; position[1] = bPos; error = err;
    }
    unsigned addContact(Contact* contact, unsigned limit) const override {
        Vector3 aPosWorld = body[0]->getPointInWorldSpace(position[0]);
        Vector3 bPosWorld = body[1]->getPointInWorldSpace(position[1]);
        Vector3 aToB = bPosWorld - aPosWorld;
        Vector3 normal = aToB;
        real length = normal.magnitude();
        if (length <= error) return 0;
        normal *= ((real)1) / length;
        contact->body[0] = body[0]; contact->body[1] = body[1];
        contact->contactNormal = normal;
        contact->contactPoint = (aPosWorld + bPosWorld) * (real)0.5;
        contact->penetration = length - error;
        contact->friction = 1; contact->restitution = 0;
        return 1;
    }
};

} // namespace phys
