// Strand-based hair (Verlet + distance/bend constraints). Each strand is a chain
// of particles rooted on the scalp; segment constraints keep length, bend
// constraints give stiffness/curl, and the strand collides with the head sphere.
// Gravity and wind drive the motion. Built on the engine's Vector3.
#pragma once
#include "core.h"
#include <vector>

namespace phys {

class Hair {
public:
    struct Strand {
        std::vector<Vector3> pos, prev;
        std::vector<real> segLen;              // per-segment rest lengths
        std::vector<real> bendRest;            // rest distance i..i+2
        real widthRoot, widthTip;              // for rendering
        real shade;                            // per-strand colour jitter
    };
    std::vector<Strand> strands;

    Vector3 gravity = Vector3(0, -9.81, 0);
    real damping = (real)0.97;
    int iterations = 16;
    real bendStiffness = (real)0.85;           // 0 = floppy, 1 = stiff
    Vector3 headCentre = Vector3(0, 0, 0);
    real headRadius = 1;
    real collisionMargin = (real)0.01;
    real turbulence = 0;                        // per-strand fluttering strength
    real time = 0;                              // internal clock for the turbulence field

    // Add a strand growing from `root` along `dir`, n segments of length `seg`.
    // A little downward+outward droop is baked into the rest shape.
    void addStrand(const Vector3& root, const Vector3& dir, int n, real seg,
                   real wRoot, real wTip, real shade) {
        Strand s; s.widthRoot = wRoot; s.widthTip = wTip; s.shade = shade;
        Vector3 p = root, d = dir.unit();
        s.pos.push_back(p); s.prev.push_back(p);
        Vector3 growth = d;
        for (int i = 0; i < n; i++) {
            growth = (growth + Vector3(0, -0.35, 0)).unit();     // droop as it grows
            p += growth * seg;
            s.pos.push_back(p); s.prev.push_back(p);
            s.segLen.push_back(seg);
        }
        for (int i = 0; i + 2 < (int)s.pos.size(); i++) s.bendRest.push_back((s.pos[i] - s.pos[i + 2]).magnitude());
        strands.push_back(std::move(s));
    }

    void step(real dt, const Vector3& wind = Vector3()) {
        time += dt;
        for (auto& s : strands) {
            for (int i = 1; i < (int)s.pos.size(); i++) {         // root (0) stays fixed
                Vector3 w = wind;
                if (turbulence > 0) {                             // spatially/temporally varying gust field
                    const Vector3& p = s.pos[i];
                    real ph = p.x * 3.1 + p.y * 2.3 + p.z * 2.7;
                    w += Vector3(real_sin(ph + time * 7.0), real_sin(ph * 1.3 + time * 5.0 + 1.7),
                                 real_cos(ph * 0.7 + time * 6.0 + 3.1)) * turbulence;
                }
                Vector3 acc = gravity + w;
                Vector3 temp = s.pos[i];
                s.pos[i] += (s.pos[i] - s.prev[i]) * damping + acc * (dt * dt);
                s.prev[i] = temp;
            }
            for (int it = 0; it < iterations; it++) {
                for (int i = 0; i < (int)s.segLen.size(); i++) {  // segment length
                    Vector3 delta = s.pos[i + 1] - s.pos[i];
                    real len = delta.magnitude(); if (len < 1e-9) continue;
                    real diff = (len - s.segLen[i]) / len;
                    if (i == 0) s.pos[i + 1] -= delta * diff;      // root fixed
                    else { s.pos[i] += delta * (diff * (real)0.5); s.pos[i + 1] -= delta * (diff * (real)0.5); }
                }
                for (int i = 0; i < (int)s.bendRest.size(); i++) { // bending stiffness
                    Vector3 delta = s.pos[i + 2] - s.pos[i];
                    real len = delta.magnitude(); if (len < 1e-9) continue;
                    real diff = (len - s.bendRest[i]) / len * bendStiffness;
                    if (i == 0) s.pos[i + 2] -= delta * diff;
                    else { s.pos[i] += delta * (diff * (real)0.5); s.pos[i + 2] -= delta * (diff * (real)0.5); }
                }
                for (int i = 1; i < (int)s.pos.size(); i++) {      // head collision
                    Vector3 d = s.pos[i] - headCentre; real len = d.magnitude(); real r = headRadius + collisionMargin;
                    if (len < r && len > 1e-9) s.pos[i] = headCentre + d * (r / len);
                }
            }
        }
    }
};

} // namespace phys
