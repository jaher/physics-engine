// Destructible objects. A solid is pre-split into box-shaped fragments that are
// "welded" (held asleep, moving as one) until shatter() bursts them apart with a
// radial + directional impulse; from then on they are ordinary rigid bodies that
// collide and settle. Grid fracture with jittered split planes gives irregular
// shards (chunky for stone, long splinters for wood).
#pragma once
#include "core.h"
#include "body.h"
#include <vector>
#include <algorithm>

namespace phys {

struct FragmentDesc { Vector3 centre; Vector3 halfSize; };

// Split the box (centre, halfSize) into nx·ny·nz fragments. `jitter` in [0,~0.6]
// randomly offsets the interior split planes, so fragments are irregular.
inline std::vector<FragmentDesc> fractureBoxGrid(const Vector3& centre, const Vector3& halfSize,
                                                 int nx, int ny, int nz, real jitter, unsigned seed) {
    unsigned st = seed ? seed : 1u;
    auto rr = [&]() { st = st * 1103515245u + 12345u; return (real)(((st >> 16) & 0x7fff) / 32767.0); };
    auto splits = [&](int n, real half) {
        std::vector<real> b(n + 1); b[0] = -half; b[n] = half;
        for (int i = 1; i < n; i++) { real base = -half + 2 * half * i / n; b[i] = base + (rr() * 2 - 1) * jitter * (2 * half / n); }
        std::sort(b.begin(), b.end()); return b;
    };
    auto bx = splits(nx, halfSize.x), by = splits(ny, halfSize.y), bz = splits(nz, halfSize.z);
    std::vector<FragmentDesc> out;
    for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) for (int k = 0; k < nz; k++) {
        Vector3 lo(bx[i], by[j], bz[k]), hi(bx[i + 1], by[j + 1], bz[k + 1]);
        Vector3 hh = (hi - lo) * (real)0.5;
        if (hh.x < 1e-4 || hh.y < 1e-4 || hh.z < 1e-4) continue;
        out.push_back({centre + (lo + hi) * (real)0.5, hh});
    }
    return out;
}

// Holds the fragment bodies (owned elsewhere) and bursts them on demand.
struct Destructible {
    std::vector<RigidBody*> fragments;
    bool shattered = false;

    // impact: where the blow lands; forward: momentum imparted to every fragment;
    // burst: radial explosion strength (falls off with distance); spin: tumble.
    void shatter(const Vector3& impact, const Vector3& forward, real burst, real spin, unsigned seed = 7) {
        shattered = true; unsigned st = seed ? seed : 1u;
        auto rr = [&]() { st = st * 1103515245u + 12345u; return (real)(((st >> 16) & 0x7fff) / 32767.0); };
        for (auto* f : fragments) {
            Vector3 d = f->getPosition() - impact; real dist = d.magnitude() + (real)0.15;
            Vector3 dir = d * (((real)1) / dist);
            Vector3 vel = dir * (burst / dist) + forward
                        + Vector3(rr() - (real)0.5, rr() - (real)0.5, rr() - (real)0.5) * (real)0.8;
            f->setAwake(true);
            f->setVelocity(vel);
            f->setRotation(Vector3(rr() - (real)0.5, rr() - (real)0.5, rr() - (real)0.5) * spin);
        }
    }

    // Detonate an explosive charge at `centre`. Every fragment within `radius` is
    // thrown radially outward at a speed that falls off with distance (blast
    // overpressure ~ 1/r in the near field), plus an upward plume bias and a random
    // tumble; fragments beyond the radius are left undisturbed. `strength` sets the
    // blast velocity scale, `upBias` the fraction added as a rising plume.
    void detonate(const Vector3& centre, real strength, real upBias, real radius, real spin, unsigned seed = 7) {
        shattered = true; unsigned st = seed ? seed : 1u;
        auto rr = [&]() { st = st * 1103515245u + 12345u; return (real)(((st >> 16) & 0x7fff) / 32767.0); };
        const real core = (real)0.25;                       // near-field cap so v stays finite at the charge
        for (auto* f : fragments) {
            Vector3 d = f->getPosition() - centre; real dist = d.magnitude();
            if (dist > radius) continue;                    // outside the blast — untouched
            Vector3 dir = dist > (real)1e-6 ? d * (((real)1) / dist) : Vector3(0, 1, 0);
            real speed = strength / (dist + core);          // overpressure falloff
            Vector3 vel = dir * speed
                        + Vector3(0, 1, 0) * (upBias * speed)
                        + Vector3(rr() - (real)0.5, rr() - (real)0.5, rr() - (real)0.5) * (strength * (real)0.12);
            f->setAwake(true);
            f->setVelocity(vel);
            f->setRotation(Vector3(rr() - (real)0.5, rr() - (real)0.5, rr() - (real)0.5) * spin);
        }
    }
};

} // namespace phys
