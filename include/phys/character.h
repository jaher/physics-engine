// Kinematic character controller (PhysX/Unity-style): a capsule moved with
// collide-and-slide against planes, boxes, spheres and heightfield terrain, with
// gravity, jumping and a grounded flag. Not force-driven — positional, like the
// controllers games actually use.
#pragma once
#include "shapes.h"
#include "terrain.h"
#include <vector>

namespace phys {

class CharacterController {
public:
    Vector3 position;                                      // capsule centre
    real radius = (real)0.35, halfHeight = (real)0.55;     // capsule spine half-length
    Vector3 velocity;
    real gravity = (real)-9.81;
    bool grounded = false;
    real maxSlopeCos = (real)0.5;                          // steeper than ~60° is a wall

    // static world
    struct SPlane { Vector3 n; real offset; };
    std::vector<SPlane> planes;
    std::vector<CollisionBox*> boxes;
    struct SSphere { Vector3 c; real r; };
    std::vector<SSphere> spheres;
    const HeightField* terrain = nullptr;

    void jump(real speed) { if (grounded) { velocity.y = speed; grounded = false; } }

    // one controller step: desiredVel is the wish velocity in the ground plane
    void move(const Vector3& desiredVel, real dt) {
        velocity.x = desiredVel.x; velocity.z = desiredVel.z;
        if (!grounded) velocity.y += gravity * dt;
        else if (velocity.y < 0) velocity.y = 0;
        position += velocity * dt;
        grounded = false;
        for (int pass = 0; pass < 4; pass++) resolveOnce();
    }
    real skin = (real)0.02;                                // contact tolerance (prevents grounded flicker)
private:
    void push(const Vector3& n, real pen) {
        if (pen <= -skin) return;
        if (n.y > maxSlopeCos) grounded = true;            // within the skin counts as standing
        if (pen <= 0) return;
        position += n * pen;
        real vn = velocity * n; if (vn < 0) velocity -= n * vn;   // slide: kill inward velocity
    }
    void spineEndpoints(Vector3& a, Vector3& b) const {
        a = position + Vector3(0, halfHeight, 0); b = position - Vector3(0, halfHeight, 0);
    }
    void resolveOnce() {
        Vector3 a, b; spineEndpoints(a, b);
        for (auto& pl : planes) {                          // both capsule end-spheres
            for (const Vector3& e : {a, b}) { real d = pl.n * e - pl.offset;
                if (d < radius + skin) push(pl.n, radius - d); spineEndpoints(a, b); }
        }
        for (auto& sp : spheres) {
            Vector3 p = closestOnSegment(a, b, sp.c);
            Vector3 d = p - sp.c; real len = d.magnitude(); real R = radius + sp.r;
            if (len < R + skin && len > 1e-9) { push(d * (((real)1) / len), R - len); spineEndpoints(a, b); }
        }
        for (auto* bx : boxes) {
            bx->calculateInternals();
            for (const Vector3& e : {a, b}) {
                Vector3 rel = bx->getTransform().transformInverse(e);
                Vector3 cl = rel;
                for (int i = 0; i < 3; i++) { if (cl[i] > bx->halfSize[i]) cl[i] = bx->halfSize[i];
                    if (cl[i] < -bx->halfSize[i]) cl[i] = -bx->halfSize[i]; }
                Vector3 cw = bx->getTransform().transform(cl);
                Vector3 d = e - cw; real len = d.magnitude();
                if (len < radius + skin && len > 1e-9) { push(d * (((real)1) / len), radius - len); spineEndpoints(a, b); }
            }
        }
        if (terrain) {
            Vector3 foot = position - Vector3(0, halfHeight, 0);
            real ground = terrain->sample(foot.x, foot.z);
            real d = foot.y - ground;
            if (d < radius + skin) { Vector3 n = terrain->normal(foot.x, foot.z); push(n, (radius - d) * n.y); }
        }
    }
};

} // namespace phys
