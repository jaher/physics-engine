// Coarse (broad-phase) collision detection (Millington Part IV ch.12): a
// bounding-sphere hierarchy (BVH) that quickly finds pairs of bodies whose
// bounding volumes overlap, so the expensive fine tests run only on candidates.
#pragma once
#include "body.h"
#include <vector>

namespace phys {

struct BoundingSphere {
    Vector3 centre;
    real radius;
    BoundingSphere() : radius(0) {}
    BoundingSphere(const Vector3& centre, real radius) : centre(centre), radius(radius) {}
    // smallest sphere enclosing two others
    BoundingSphere(const BoundingSphere& one, const BoundingSphere& two) {
        Vector3 centreOffset = two.centre - one.centre;
        real distanceSq = centreOffset.squareMagnitude();
        real radiusDiff = two.radius - one.radius;
        if (radiusDiff * radiusDiff >= distanceSq) {          // one contains the other
            if (one.radius > two.radius) { centre = one.centre; radius = one.radius; }
            else { centre = two.centre; radius = two.radius; }
        } else {
            real distance = real_sqrt(distanceSq);
            radius = (distance + one.radius + two.radius) * (real)0.5;
            centre = one.centre;
            if (distance > 0) centre += centreOffset * ((radius - one.radius) / distance);
        }
    }
    bool overlaps(const BoundingSphere* other) const {
        real distanceSq = (centre - other->centre).squareMagnitude();
        return distanceSq < (radius + other->radius) * (radius + other->radius);
    }
    real getGrowth(const BoundingSphere& other) const {
        BoundingSphere newSphere(*this, other);
        return newSphere.radius * newSphere.radius - radius * radius;
    }
    real getSize() const { return ((real)1.333333) * real_pi * radius * radius * radius; }
};

struct PotentialContact { RigidBody* body[2]; };

template <class BoundingVolumeClass>
class BVHNode {
public:
    BVHNode* children[2];
    BoundingVolumeClass volume;
    RigidBody* body;
    BVHNode* parent;

    BVHNode(BVHNode* parent, const BoundingVolumeClass& volume, RigidBody* body = nullptr)
        : volume(volume), body(body), parent(parent) { children[0] = children[1] = nullptr; }
    ~BVHNode() {
        if (parent) {
            BVHNode* sibling = (parent->children[0] == this) ? parent->children[1] : parent->children[0];
            parent->volume = sibling->volume; parent->body = sibling->body;
            parent->children[0] = sibling->children[0]; parent->children[1] = sibling->children[1];
            if (sibling->children[0]) sibling->children[0]->parent = parent;
            if (sibling->children[1]) sibling->children[1]->parent = parent;
            sibling->parent = nullptr; sibling->body = nullptr;
            sibling->children[0] = sibling->children[1] = nullptr;
            delete sibling;
            parent->recalculateBoundingVolume();
        }
        if (children[0]) { children[0]->parent = nullptr; delete children[0]; }
        if (children[1]) { children[1]->parent = nullptr; delete children[1]; }
    }
    bool isLeaf() const { return body != nullptr; }
    bool overlaps(const BVHNode* other) const { return volume.overlaps(&other->volume); }

    unsigned getPotentialContacts(PotentialContact* contacts, unsigned limit) const {
        if (isLeaf() || limit == 0) return 0;
        return children[0]->getPotentialContactsWith(children[1], contacts, limit);
    }
    unsigned getPotentialContactsWith(const BVHNode* other, PotentialContact* contacts, unsigned limit) const {
        if (!overlaps(other) || limit == 0) return 0;
        if (isLeaf() && other->isLeaf()) { contacts->body[0] = body; contacts->body[1] = other->body; return 1; }
        if (other->isLeaf() || (!isLeaf() && volume.getSize() >= other->volume.getSize())) {
            unsigned count = children[0]->getPotentialContactsWith(other, contacts, limit);
            if (limit > count) count += children[1]->getPotentialContactsWith(other, contacts + count, limit - count);
            return count;
        } else {
            unsigned count = getPotentialContactsWith(other->children[0], contacts, limit);
            if (limit > count) count += getPotentialContactsWith(other->children[1], contacts + count, limit - count);
            return count;
        }
    }
    void insert(RigidBody* newBody, const BoundingVolumeClass& newVolume) {
        if (isLeaf()) {
            children[0] = new BVHNode(this, volume, body);
            children[1] = new BVHNode(this, newVolume, newBody);
            body = nullptr;
            recalculateBoundingVolume();
        } else {
            if (children[0]->volume.getGrowth(newVolume) < children[1]->volume.getGrowth(newVolume))
                children[0]->insert(newBody, newVolume);
            else
                children[1]->insert(newBody, newVolume);
        }
    }
    void recalculateBoundingVolume(bool recurse = true) {
        if (isLeaf()) return;
        volume = BoundingVolumeClass(children[0]->volume, children[1]->volume);
        if (parent && recurse) parent->recalculateBoundingVolume(true);
    }
};

} // namespace phys
