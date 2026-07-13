// Compound (multi-shape) rigid bodies — a single body carrying several collision
// primitives at local offsets, with the aggregate mass, centre of mass and
// inertia tensor computed by the parallel-axis theorem. Every mainstream engine
// has this; the primitives already support a local `offset`, so this adds the
// mass-property bookkeeping and builds the offset primitives for the contact pass.
#pragma once
#include "collide_fine.h"
#include <vector>
#include <memory>

namespace phys {

class CompoundShape {
public:
    struct Child { bool box; Vector3 half; real radius; Vector3 pos; Quaternion orient; real mass; };
    std::vector<Child> children;
    real totalMass = 0; Vector3 com; Matrix3 inertia;              // aggregate mass properties (about the COM)
    std::vector<std::unique_ptr<CollisionBox>> boxes;
    std::vector<std::unique_ptr<CollisionSphere>> spheres;

    void addBox(const Vector3& half, const Vector3& pos, real mass, const Quaternion& o = Quaternion()) { children.push_back({true, half, 0, pos, o, mass}); }
    void addSphere(real radius, const Vector3& pos, real mass) { children.push_back({false, Vector3(), radius, pos, Quaternion(), mass}); }

    // Aggregate the children onto `body`: total mass, COM, inertia (each child's
    // inertia rotated into the body frame + parallel-axis shifted to the COM), and
    // build the offset collision primitives (relative to the COM = body origin).
    void finalize(RigidBody* body) {
        totalMass = 0; com = Vector3();
        for (auto& c : children) { totalMass += c.mass; com += c.pos * c.mass; }
        if (totalMass > 0) com *= ((real)1) / totalMass;
        inertia = Matrix3();                                        // zeroed
        for (auto& c : children) {
            Matrix3 Ic;
            if (c.box) { Ic.setBlockInertiaTensor(c.half, c.mass); Matrix3 R; R.setOrientation(c.orient); Ic = R * Ic * R.transpose(); }
            else { real v = (real)0.4 * c.mass * c.radius * c.radius; Ic.setInertiaTensorCoeffs(v, v, v); }
            Vector3 d = c.pos - com; real d2 = d * d;               // parallel axis: I += m(|d|²E − d⊗d)
            Matrix3 pa;
            pa.data[0] = c.mass * (d2 - d.x * d.x); pa.data[1] = -c.mass * d.x * d.y; pa.data[2] = -c.mass * d.x * d.z;
            pa.data[3] = -c.mass * d.y * d.x; pa.data[4] = c.mass * (d2 - d.y * d.y); pa.data[5] = -c.mass * d.y * d.z;
            pa.data[6] = -c.mass * d.z * d.x; pa.data[7] = -c.mass * d.z * d.y; pa.data[8] = c.mass * (d2 - d.z * d.z);
            inertia += Ic; inertia += pa;
        }
        body->setMass(totalMass); body->setInertiaTensor(inertia);
        boxes.clear(); spheres.clear();
        for (auto& c : children) {
            Matrix4 off; off.setOrientationAndPos(c.orient, c.pos - com);
            if (c.box) { auto b = std::make_unique<CollisionBox>(); b->body = body; b->halfSize = c.half; b->offset = off; boxes.push_back(std::move(b)); }
            else { auto s = std::make_unique<CollisionSphere>(); s->body = body; s->radius = c.radius; s->offset = off; spheres.push_back(std::move(s)); }
        }
    }
    // update every child primitive's world transform (call before contact generation)
    void calculateInternals() { for (auto& b : boxes) b->calculateInternals(); for (auto& s : spheres) s->calculateInternals(); }
};

} // namespace phys
