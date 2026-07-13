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
};

} // namespace phys
