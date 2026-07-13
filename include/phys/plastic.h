// Elastoplastic rod with thermal softening — iron you can heat and bend. The rod
// is a Verlet particle chain with distance + bend constraints. Each bend
// constraint stores its REST configuration; when the elastic strain exceeds a
// yield threshold the rest state creeps toward the current one (plastic flow,
// permanent set). Heat lowers the yield threshold and softens the stiffness, so
// a cold rod springs back where a red-hot one takes the bend. Temperature
// diffuses along the rod and radiates away.
#pragma once
#include "core.h"
#include <vector>

namespace phys {

class PlasticRod {
public:
    std::vector<Vector3> pos, prev;
    std::vector<unsigned char> pinned;
    std::vector<real> segRest;                 // rest length per segment
    struct Bend { int a, b; real rest; };      // multi-scale bend constraints (spans 2,4,8)
    std::vector<Bend> bends;
    std::vector<real> temp;                    // per-particle temperature, 0 = cold, 1 = forging heat

    Vector3 gravity = Vector3(0, -9.81, 0);
    real damping = (real)0.94;
    int iterations = 24;
    real stiffCold = (real)0.95, stiffHot = (real)0.08;    // bend stiffness vs temperature
    real yieldCold = (real)0.05, yieldHot = (real)0.0003;  // elastic strain limit vs temperature
    real plasticRate = (real)0.9;              // how fast rest state creeps once yielding
    real coolRate = (real)0.05;                // radiative cooling per second
    real conduct = (real)0.8;                  // conduction along the rod per second

    void build(const Vector3& start, const Vector3& dir, int n, real seg) {
        pos.clear(); prev.clear(); pinned.clear(); segRest.clear(); bends.clear(); temp.clear();
        Vector3 d = dir.unit(), p = start;
        for (int i = 0; i <= n; i++) { pos.push_back(p); prev.push_back(p); pinned.push_back(0);
            temp.push_back(0); if (i < n) segRest.push_back(seg); p += d * seg; }
        for (int span : {2, 4, 8})
            for (int i = 0; i + span < (int)pos.size(); i++)
                bends.push_back({i, i + span, (pos[i] - pos[i + span]).magnitude()});
    }
    void pin(int i) { pinned[i] = 1; }
    // pour heat into the rod around arc position `at` (particle index space)
    void heat(real at, real radius, real amount, real dt) {
        for (int i = 0; i < (int)temp.size(); i++) {
            real d = real_abs((real)i - at) / radius;
            if (d < 1) { temp[i] += amount * (1 - d) * dt; if (temp[i] > 1) temp[i] = 1; }
        }
    }
    real maxTemp() const { real m = 0; for (real t : temp) if (t > m) m = t; return m; }

    void step(real dt, const Vector3& tipForce = Vector3(), int tipIndex = -1) {
        // integrate
        for (int i = 0; i < (int)pos.size(); i++) {
            if (pinned[i]) { prev[i] = pos[i]; continue; }
            Vector3 acc = gravity;
            if (i == tipIndex) acc += tipForce;            // external tool push (per unit mass)
            Vector3 tmp = pos[i];
            pos[i] += (pos[i] - prev[i]) * damping + acc * (dt * dt);
            prev[i] = tmp;
        }
        // solve constraints
        for (int it = 0; it < iterations; it++) {
            for (int i = 0; i < (int)segRest.size(); i++) {               // inextensible segments
                Vector3 d = pos[i + 1] - pos[i]; real len = d.magnitude(); if (len < 1e-9) continue;
                real diff = (len - segRest[i]) / len;
                bool pa = pinned[i], pb = pinned[i + 1];
                if (pa && pb) continue;
                if (!pa && !pb) { pos[i] += d * (diff * (real)0.5); pos[i + 1] -= d * (diff * (real)0.5); }
                else if (pa) pos[i + 1] -= d * diff;
                else pos[i] += d * diff;
            }
            for (auto& bc : bends) {                                     // multi-scale bend toward rest
                real t = (temp[bc.a] + temp[bc.b]) / 2;
                real k = stiffCold + (stiffHot - stiffCold) * t;          // hot → soft
                Vector3 d = pos[bc.b] - pos[bc.a]; real len = d.magnitude(); if (len < 1e-9) continue;
                real diff = (len - bc.rest) / len * k;
                bool pa = pinned[bc.a], pc = pinned[bc.b];
                if (pa && pc) continue;
                if (!pa && !pc) { pos[bc.a] += d * (diff * (real)0.5); pos[bc.b] -= d * (diff * (real)0.5); }
                else if (pa) pos[bc.b] -= d * diff;
                else pos[bc.a] += d * diff;
            }
        }
        // plastic flow: where elastic strain exceeds the (temperature-dependent)
        // yield, the rest configuration creeps toward the current one.
        for (auto& bc : bends) {
            real t = (temp[bc.a] + temp[bc.b]) / 2;
            real yield = yieldCold + (yieldHot - yieldCold) * t;
            real cur = (pos[bc.a] - pos[bc.b]).magnitude();
            real strain = (cur - bc.rest) / bc.rest;
            if (real_abs(strain) > yield) {
                real target = cur - (strain > 0 ? yield : -yield) * bc.rest;
                bc.rest += (target - bc.rest) * plasticRate * (t > (real)0.15 ? 1 : (real)0.05);
            }
        }
        // thermal: conduction + radiative cooling
        std::vector<real> nt = temp;
        for (int i = 0; i < (int)temp.size(); i++) {
            real lap = 0; int n = 0;
            if (i > 0) { lap += temp[i - 1] - temp[i]; n++; }
            if (i + 1 < (int)temp.size()) { lap += temp[i + 1] - temp[i]; n++; }
            nt[i] += conduct * dt * lap / (n ? n : 1);
            nt[i] -= coolRate * dt * nt[i];
            if (nt[i] < 0) nt[i] = 0;
        }
        temp = nt;
    }
    // permanent deflection of the tip from the straight rest line (for tests)
    real permanentSetAngle() const {
        real bent = 0;
        for (auto& bc : bends) {
            if (bc.b - bc.a != 2) continue;                 // measure from span-2 constraints
            real s0 = segRest[bc.a], s1 = segRest[bc.a + 1];
            real c = (s0 * s0 + s1 * s1 - bc.rest * bc.rest) / (2 * s0 * s1);
            if (c > 1) c = 1; if (c < -1) c = -1;
            bent += real_pi - std::acos((double)c);
        }
        return bent;
    }
};

} // namespace phys
