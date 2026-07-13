// Raycast vehicle (the Bullet raycastVehicle / Unity WheelCollider / PhysX
// vehicle model): a rigid chassis with per-wheel suspension raycasts, spring +
// damper forces, engine/brake drive and lateral tyre friction. Drives on the
// ground plane or a HeightField.
#pragma once
#include "body.h"
#include "terrain.h"
#include <vector>

namespace phys {

class RaycastVehicle {
public:
    RigidBody* chassis = nullptr;
    const HeightField* terrain = nullptr;                  // optional; else ground plane y=0
    struct Wheel {
        Vector3 localAttach;                                // suspension top, chassis space
        real restLength = (real)0.55, radius = (real)0.3;
        real stiffness = 32000, damping = 3800;
        bool steerable = false, driven = false;
        // runtime
        bool grounded = false; real compression = 0, lastCompression = 0;
        Vector3 contactPoint, contactNormal;
    };
    std::vector<Wheel> wheels;
    real steer = 0;                                        // radians
    real engineForce = 0, brakeForce = 0;
    real lateralGrip = (real)5.5;                          // tyre lateral stiffness (per unit load)
    real maxGripForce = 9000;

    void addWheel(const Vector3& attach, bool steerable, bool driven) {
        Wheel w; w.localAttach = attach; w.steerable = steerable; w.driven = driven; wheels.push_back(w);
    }
    real groundHeight(real x, real z) const { return terrain ? terrain->sample(x, z) : 0; }
    Vector3 groundNormal(real x, real z) const { return terrain ? terrain->normal(x, z) : Vector3(0, 1, 0); }

    void update(real dt) {
        Vector3 up = chassis->getDirectionInWorldSpace(Vector3(0, 1, 0));
        Vector3 fwd = chassis->getDirectionInWorldSpace(Vector3(0, 0, 1));
        for (auto& w : wheels) {
            Vector3 top = chassis->getPointInWorldSpace(w.localAttach);
            // cast straight down in world space (robust for mild chassis tilt)
            real ground = groundHeight(top.x, top.z);
            real drop = top.y - ground;                     // distance to terrain below the attach
            real maxLen = w.restLength + w.radius;
            w.lastCompression = w.compression;
            if (drop < maxLen) {
                w.grounded = true;
                w.compression = maxLen - drop;
                w.contactNormal = groundNormal(top.x, top.z);
                w.contactPoint = Vector3(top.x, ground, top.z);
                real compVel = (w.compression - w.lastCompression) / dt;
                real force = w.stiffness * w.compression + w.damping * compVel;
                if (force < 0) force = 0;
                chassis->addForceAtPoint(w.contactNormal * force, top);

                // wheel-plane directions
                Vector3 wf = fwd;
                if (w.steerable) {                          // rotate forward by the steer angle about up
                    real c = real_cos(steer), s = real_sin(steer);
                    Vector3 right = fwd % up;
                    wf = (fwd * c + right * -s).unit();
                }
                Vector3 n = w.contactNormal;
                wf = (wf - n * (wf * n)).unit();            // project into the contact plane
                Vector3 side = (n % wf).unit();

                // velocity of the contact point
                Vector3 rel = w.contactPoint - chassis->getPosition();
                Vector3 vel = chassis->getVelocity() + (chassis->getRotation() % rel);
                real latVel = vel * side, fwdVel = vel * wf;

                // lateral tyre friction (proportional to load, clamped)
                real lat = -latVel * lateralGrip * force * dt * 60 / 60;
                real cap = maxGripForce; if (lat > cap) lat = cap; if (lat < -cap) lat = -cap;
                chassis->addForceAtPoint(side * lat, w.contactPoint);

                // drive + brake
                if (w.driven && engineForce != 0) chassis->addForceAtPoint(wf * engineForce, w.contactPoint);
                if (brakeForce > 0) { real br = -fwdVel * brakeForce; if (br > brakeForce) br = brakeForce;
                    if (br < -brakeForce) br = -brakeForce; chassis->addForceAtPoint(wf * br, w.contactPoint); }
            } else { w.grounded = false; w.compression = 0; }
        }
    }
    // world transform of a wheel's hub (for rendering)
    Vector3 wheelCentre(const Wheel& w) const {
        Vector3 top = chassis->getPointInWorldSpace(w.localAttach);
        real ground = groundHeight(top.x, top.z);
        real drop = top.y - ground;
        real len = w.grounded ? (drop - w.radius) : w.restLength;
        return top - Vector3(0, len, 0);
    }
};

} // namespace phys
