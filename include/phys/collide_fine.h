// Fine (narrow-phase) collision detection (Millington Part IV ch.13). Collision
// primitives (sphere, box, half-space plane), quick boolean intersection tests,
// and a detector that generates full Contacts: sphere-sphere, sphere-plane,
// box-plane, box-sphere and box-box (separating-axis test with face and
// edge-edge contact generation).
#pragma once
#include "contacts.h"
#include <cassert>

namespace phys {

class CollisionPrimitive {
public:
    RigidBody* body = nullptr;
    Matrix4 offset;                        // primitive offset from the body origin
    Matrix4 transform;                     // derived world transform (body * offset)
    void calculateInternals() { transform = body->getTransform() * offset; }
    Vector3 getAxis(unsigned i) const { return transform.getAxisVector(i); }
    const Matrix4& getTransform() const { return transform; }
};

class CollisionSphere : public CollisionPrimitive { public: real radius = 1; };
class CollisionBox : public CollisionPrimitive { public: Vector3 halfSize; };

// An immovable half-space (points p with p·direction ≤ offset are "inside").
class CollisionPlane { public: Vector3 direction; real offset = 0; };

struct CollisionData {
    Contact* contactArray = nullptr;       // start of the writable array
    Contact* contacts = nullptr;           // current write head
    int contactsLeft = 0;
    unsigned contactCount = 0;
    real friction = 0, restitution = 0, tolerance = 0;
    bool hasMoreContacts() const { return contactsLeft > 0; }
    void reset(unsigned maxContacts) { contactsLeft = maxContacts; contactCount = 0; contacts = contactArray; }
    void addContacts(unsigned count) { contactsLeft -= count; contactCount += count; contacts += count; }
};

inline real transformToAxis(const CollisionBox& box, const Vector3& axis) {
    return box.halfSize.x * real_abs(axis * box.getAxis(0))
         + box.halfSize.y * real_abs(axis * box.getAxis(1))
         + box.halfSize.z * real_abs(axis * box.getAxis(2));
}

class IntersectionTests {
public:
    static bool sphereAndHalfSpace(const CollisionSphere& sphere, const CollisionPlane& plane) {
        real ballDistance = plane.direction * sphere.getAxis(3) - sphere.radius;
        return ballDistance <= plane.offset;
    }
    static bool sphereAndSphere(const CollisionSphere& one, const CollisionSphere& two) {
        Vector3 midline = one.getAxis(3) - two.getAxis(3);
        return midline.squareMagnitude() < (one.radius + two.radius) * (one.radius + two.radius);
    }
    static bool boxAndHalfSpace(const CollisionBox& box, const CollisionPlane& plane) {
        real projectedRadius = transformToAxis(box, plane.direction);
        real boxDistance = plane.direction * box.getAxis(3) - projectedRadius;
        return boxDistance <= plane.offset;
    }
};

class CollisionDetector {
    static real penetrationOnAxis(const CollisionBox& one, const CollisionBox& two, const Vector3& axis, const Vector3& toCentre) {
        return transformToAxis(one, axis) + transformToAxis(two, axis) - real_abs(toCentre * axis);
    }
    static bool tryAxis(const CollisionBox& one, const CollisionBox& two, Vector3 axis, const Vector3& toCentre,
                        unsigned index, real& smallestPenetration, unsigned& smallestCase) {
        if (axis.squareMagnitude() < (real)0.0001) return true;   // near-parallel edges → skip
        axis.normalise();
        real penetration = penetrationOnAxis(one, two, axis, toCentre);
        if (penetration < 0) return false;
        if (penetration < smallestPenetration) { smallestPenetration = penetration; smallestCase = index; }
        return true;
    }
    static void fillPointFaceBoxBox(const CollisionBox& one, const CollisionBox& two, const Vector3& toCentre,
                                    CollisionData* data, unsigned best, real pen) {
        Contact* contact = data->contacts;
        Vector3 normal = one.getAxis(best);
        if (one.getAxis(best) * toCentre > 0) normal *= -1;
        Vector3 vertex = two.halfSize;
        if (two.getAxis(0) * normal < 0) vertex.x = -vertex.x;
        if (two.getAxis(1) * normal < 0) vertex.y = -vertex.y;
        if (two.getAxis(2) * normal < 0) vertex.z = -vertex.z;
        contact->contactNormal = normal;
        contact->penetration = pen;
        contact->contactPoint = two.getTransform() * vertex;
        contact->setBodyData(one.body, two.body, data->friction, data->restitution);
    }
    static Vector3 contactPoint(const Vector3& pOne, const Vector3& dOne, real oneSize,
                                const Vector3& pTwo, const Vector3& dTwo, real twoSize, bool useOne) {
        real smOne = dOne.squareMagnitude(), smTwo = dTwo.squareMagnitude();
        real dpOneTwo = dTwo * dOne;
        Vector3 toSt = pOne - pTwo;
        real dpStaOne = dOne * toSt, dpStaTwo = dTwo * toSt;
        real denom = smOne * smTwo - dpOneTwo * dpOneTwo;
        if (real_abs(denom) < (real)0.0001) return useOne ? pOne : pTwo;
        real mua = (dpOneTwo * dpStaTwo - smTwo * dpStaOne) / denom;
        real mub = (smOne * dpStaTwo - dpOneTwo * dpStaOne) / denom;
        if (mua > oneSize || mua < -oneSize || mub > twoSize || mub < -twoSize) return useOne ? pOne : pTwo;
        Vector3 cOne = pOne + dOne * mua, cTwo = pTwo + dTwo * mub;
        return cOne * (real)0.5 + cTwo * (real)0.5;
    }
public:
    static unsigned sphereAndHalfSpace(const CollisionSphere& sphere, const CollisionPlane& plane, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 position = sphere.getAxis(3);
        real ballDistance = plane.direction * position - sphere.radius - plane.offset;
        if (ballDistance >= 0) return 0;
        Contact* c = data->contacts;
        c->contactNormal = plane.direction;
        c->penetration = -ballDistance;
        c->contactPoint = position - plane.direction * (ballDistance + sphere.radius);
        c->setBodyData(sphere.body, nullptr, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    static unsigned sphereAndTruePlane(const CollisionSphere& sphere, const CollisionPlane& plane, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 position = sphere.getAxis(3);
        real centreDistance = plane.direction * position - plane.offset;
        if (centreDistance * centreDistance > sphere.radius * sphere.radius) return 0;
        Vector3 normal = plane.direction; real penetration = -centreDistance;
        if (centreDistance < 0) { normal *= -1; penetration = -penetration; }
        penetration += sphere.radius;
        Contact* c = data->contacts;
        c->contactNormal = normal;
        c->penetration = penetration;
        c->contactPoint = position - plane.direction * centreDistance;
        c->setBodyData(sphere.body, nullptr, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    static unsigned sphereAndSphere(const CollisionSphere& one, const CollisionSphere& two, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 pos1 = one.getAxis(3), pos2 = two.getAxis(3);
        Vector3 midline = pos1 - pos2; real size = midline.magnitude();
        if (size <= 0 || size >= one.radius + two.radius) return 0;
        Vector3 normal = midline * (((real)1) / size);        // points from two toward one
        real penetration = one.radius + two.radius - size;
        Contact* c = data->contacts;
        c->contactNormal = normal;
        // geometric midpoint of the overlap region on the centre line
        c->contactPoint = pos1 - normal * (one.radius - penetration * (real)0.5);
        c->penetration = penetration;
        c->setBodyData(one.body, two.body, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    static unsigned boxAndHalfSpace(const CollisionBox& box, const CollisionPlane& plane, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        if (!IntersectionTests::boxAndHalfSpace(box, plane)) return 0;
        static real mults[8][3] = {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
        Contact* contact = data->contacts; unsigned contactsUsed = 0;
        for (unsigned i = 0; i < 8; i++) {
            Vector3 vertexPos(mults[i][0], mults[i][1], mults[i][2]);
            vertexPos.componentProductUpdate(box.halfSize);
            vertexPos = box.getTransform().transform(vertexPos);
            real vertexDistance = vertexPos * plane.direction;
            if (vertexDistance <= plane.offset) {
                contact->contactPoint = plane.direction * (vertexDistance - plane.offset) + vertexPos;
                contact->contactNormal = plane.direction;
                contact->penetration = plane.offset - vertexDistance;
                contact->setBodyData(box.body, nullptr, data->friction, data->restitution);
                contact++; contactsUsed++;
                if ((int)contactsUsed == data->contactsLeft) break;
            }
        }
        data->addContacts(contactsUsed); return contactsUsed;
    }
    static unsigned boxAndSphere(const CollisionBox& box, const CollisionSphere& sphere, CollisionData* data) {
        if (data->contactsLeft <= 0) return 0;
        Vector3 centre = sphere.getAxis(3);
        Vector3 relCentre = box.getTransform().transformInverse(centre);
        if (real_abs(relCentre.x) - sphere.radius > box.halfSize.x ||
            real_abs(relCentre.y) - sphere.radius > box.halfSize.y ||
            real_abs(relCentre.z) - sphere.radius > box.halfSize.z) return 0;
        Vector3 closestPt; real dist;
        dist = relCentre.x; if (dist > box.halfSize.x) dist = box.halfSize.x; if (dist < -box.halfSize.x) dist = -box.halfSize.x; closestPt.x = dist;
        dist = relCentre.y; if (dist > box.halfSize.y) dist = box.halfSize.y; if (dist < -box.halfSize.y) dist = -box.halfSize.y; closestPt.y = dist;
        dist = relCentre.z; if (dist > box.halfSize.z) dist = box.halfSize.z; if (dist < -box.halfSize.z) dist = -box.halfSize.z; closestPt.z = dist;
        dist = (closestPt - relCentre).squareMagnitude();
        if (dist > sphere.radius * sphere.radius) return 0;
        Vector3 closestPtWorld = box.getTransform().transform(closestPt);
        Contact* c = data->contacts;
        c->contactNormal = closestPtWorld - centre; c->contactNormal.normalise();
        c->contactPoint = closestPtWorld;
        c->penetration = sphere.radius - real_sqrt(dist);
        c->setBodyData(box.body, sphere.body, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    static unsigned boxAndBox(const CollisionBox& one, const CollisionBox& two, CollisionData* data) {
        Vector3 toCentre = two.getAxis(3) - one.getAxis(3);
        real pen = REAL_MAX; unsigned best = 0xffffff;
#define TRY(axis, index) if (!tryAxis(one, two, (axis), toCentre, (index), pen, best)) return 0;
        TRY(one.getAxis(0), 0); TRY(one.getAxis(1), 1); TRY(one.getAxis(2), 2);
        TRY(two.getAxis(0), 3); TRY(two.getAxis(1), 4); TRY(two.getAxis(2), 5);
        unsigned bestSingleAxis = best;
        TRY(one.getAxis(0) % two.getAxis(0), 6);  TRY(one.getAxis(0) % two.getAxis(1), 7);  TRY(one.getAxis(0) % two.getAxis(2), 8);
        TRY(one.getAxis(1) % two.getAxis(0), 9);  TRY(one.getAxis(1) % two.getAxis(1), 10); TRY(one.getAxis(1) % two.getAxis(2), 11);
        TRY(one.getAxis(2) % two.getAxis(0), 12); TRY(one.getAxis(2) % two.getAxis(1), 13); TRY(one.getAxis(2) % two.getAxis(2), 14);
#undef TRY
        assert(best != 0xffffff);
        if (best < 3) { fillPointFaceBoxBox(one, two, toCentre, data, best, pen); data->addContacts(1); return 1; }
        if (best < 6) { fillPointFaceBoxBox(two, one, toCentre * -1, data, best - 3, pen); data->addContacts(1); return 1; }
        // edge-edge
        best -= 6;
        unsigned oneAxisIndex = best / 3, twoAxisIndex = best % 3;
        Vector3 oneAxis = one.getAxis(oneAxisIndex), twoAxis = two.getAxis(twoAxisIndex);
        Vector3 axis = oneAxis % twoAxis; axis.normalise();
        if (axis * toCentre > 0) axis = axis * -1;
        Vector3 ptOnOneEdge = one.halfSize, ptOnTwoEdge = two.halfSize;
        for (unsigned i = 0; i < 3; i++) {
            if (i == oneAxisIndex) ptOnOneEdge[i] = 0;
            else if (one.getAxis(i) * axis > 0) ptOnOneEdge[i] = -ptOnOneEdge[i];
            if (i == twoAxisIndex) ptOnTwoEdge[i] = 0;
            else if (two.getAxis(i) * axis < 0) ptOnTwoEdge[i] = -ptOnTwoEdge[i];
        }
        ptOnOneEdge = one.getTransform() * ptOnOneEdge;
        ptOnTwoEdge = two.getTransform() * ptOnTwoEdge;
        Vector3 vertex = contactPoint(ptOnOneEdge, oneAxis, one.halfSize[oneAxisIndex],
                                      ptOnTwoEdge, twoAxis, two.halfSize[twoAxisIndex], bestSingleAxis > 2);
        Contact* contact = data->contacts;
        contact->penetration = pen;
        contact->contactNormal = axis;
        contact->contactPoint = vertex;
        contact->setBodyData(one.body, two.body, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
};

} // namespace phys
