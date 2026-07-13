// A falling, burning oak leaf: tumbling-plate aerodynamics coupled to a small
// combustion cellular automaton.
//
// Motion — the leaf is a thin plate. Anisotropic quadratic drag (strong across
// the face, weak edge-on), a reactive lift perpendicular to the airflow, and an
// aerodynamic torque that turns the leaf broadside to the flow but *overshoots*
// (underdamped) together reproduce the characteristic flutter/tumble/spiral of a
// real falling leaf. Orientation is carried as an orthonormal frame (e0,e1 in the
// leaf plane, e2 the normal) integrated by the angular velocity.
//
// Burning — cells inside an oak-leaf silhouette carry fuel/heat/char. Once lit,
// heat diffuses, consumes fuel (opening holes), blackens the cell to char, and
// curls it out of plane — the same fire model as the burning paper, per leaf.
#pragma once
#include "core.h"
#include <vector>
#include <cmath>

namespace phys {

class FallingLeaf {
public:
    int GW = 12, GH = 18;                 // grid: width × length
    real cw = 0, ch = 0;                  // physical cell size
    Vector3 pos, vel;                     // centre-of-mass state
    Vector3 e0 = Vector3(1, 0, 0), e1 = Vector3(0, 0, 1), e2 = Vector3(0, 1, 0);  // leaf frame (e2 = normal)
    Vector3 omega;                        // angular velocity (world)
    real mass = 6e-4, inertia = 2e-6, leafArea = 0;

    std::vector<unsigned char> inMask;    // oak silhouette
    std::vector<float> fuel, heat, charAmt, curl, vein, resist;   // resist → ragged, undulating burn front
    Vector3 baseColor = Vector3(0.55, 0.30, 0.12);
    int texLayer = 0;                     // which leaf texture (for the textured demo)

    // aero tuning
    real kNormal = 5.0, kTang = 0.30, kLift = 2.8, kAlign = 0.005, kAngDamp = 8e-6, airDamp = 0.999;
    // fire tuning
    real ignite = 0.5, burnRate = 0.55, spread = 1.7, cool = 0.7, curlAmt = 0.5;
    double time = 0; real igniteAt = 1e9, igniteY = -1e9; bool ignited = false, lit0 = false;

    int idx(int i, int j) const { return j * GW + i; }
    bool alive(int i, int j) const { return inMask[idx(i, j)] && fuel[idx(i, j)] > 0.04f; }
    int aliveCount = 0, maskCount = 0;
    real aliveFrac() const { return maskCount ? (real)aliveCount / maskCount : 0; }

    // oak-leaf half-width profile (v: 0 stem → 1 tip), normalised to [0,0.5]
    static real oakHalfWidth(real v) {
        real env = std::pow(std::sin(3.14159265f * std::min<real>(std::max<real>(v, 0), 1)), (real)0.55);
        real lobes = (real)0.34 + (real)0.66 * std::pow(std::fabs(std::cos(3.14159265 * 5.0 * v + 0.35)), (real)0.7);
        real hw = (real)0.47 * env * lobes;
        if (v < 0.10) hw = std::min(hw, (real)0.05 + v * (real)0.3);   // thin petiole
        return hw;
    }

    void build(real lengthM, real widthM, unsigned seed) {
        cw = widthM / GW; ch = lengthM / GH;
        int n = GW * GH;
        inMask.assign(n, 0); fuel.assign(n, 0); heat.assign(n, 0); charAmt.assign(n, 0); curl.assign(n, 0); vein.assign(n, 0); resist.assign(n, 1.0f);
        auto rr = [&]() { seed = seed * 1103515245u + 12345u; return (float)(((seed >> 16) & 0x7fff) / 32767.0f); };
        maskCount = 0;
        for (int j = 0; j < GH; j++) for (int i = 0; i < GW; i++) {
            real u = (real)i / (GW - 1), v = (real)j / (GH - 1);
            real hw = oakHalfWidth(v);
            if (std::fabs(u - 0.5) <= hw) {
                inMask[idx(i, j)] = 1; fuel[idx(i, j)] = 0.85f + 0.15f * rr(); resist[idx(i, j)] = 0.5f + 1.1f * rr(); maskCount++;
                // veins: midrib + a few lateral veins per lobe
                real midrib = 1.0 - std::min<real>(1.0, std::fabs(u - 0.5) / 0.04);
                real lat = std::pow(std::fabs(std::sin(3.14159265 * 6.0 * v)), 8.0) * (1.0 - std::fabs(u - 0.5) / std::max<real>(hw, 1e-3));
                vein[idx(i, j)] = (float)std::max<real>(midrib, lat * 0.7);
            }
        }
        leafArea = (real)maskCount * cw * ch;
        aliveCount = maskCount;
        omega = Vector3((rr() - 0.5) * 3, (rr() - 0.5) * 3, (rr() - 0.5) * 3);
    }

    // Build the leaf silhouette by sampling the alpha channel of a real leaf
    // texture (row 0 = tip, last row = stem), instead of the procedural oak mask.
    void buildTextured(int GW_, int GH_, real lengthM, real widthM,
                       const unsigned char* alpha, int tw, int th, unsigned seed) {
        GW = GW_; GH = GH_; cw = widthM / GW; ch = lengthM / GH;
        int n = GW * GH;
        inMask.assign(n, 0); fuel.assign(n, 0); heat.assign(n, 0); charAmt.assign(n, 0); curl.assign(n, 0); vein.assign(n, 0); resist.assign(n, 1.0f);
        auto rr = [&]() { seed = seed * 1103515245u + 12345u; return (float)(((seed >> 16) & 0x7fff) / 32767.0f); };
        maskCount = 0;
        for (int j = 0; j < GH; j++) for (int i = 0; i < GW; i++) {
            real u = (real)i / (GW - 1), v = (real)j / (GH - 1);
            int tx = (int)(u * (tw - 1)), ty = (int)((1.0 - v) * (th - 1));   // stem (v=0) → bottom of the image
            if (alpha[ty * tw + tx] > 128) { inMask[idx(i, j)] = 1; fuel[idx(i, j)] = 0.85f + 0.15f * rr(); resist[idx(i, j)] = 0.5f + 1.1f * rr(); maskCount++; }
        }
        leafArea = (real)maskCount * cw * ch; aliveCount = maskCount;
        omega = Vector3((rr() - 0.5) * 3, (rr() - 0.5) * 3, (rr() - 0.5) * 3);
    }

    void igniteEdge() {
        // light a random in-mask boundary cell (an edge catches first)
        int best = -1; unsigned s = (unsigned)(pos.x * 9871 + pos.z * 1237) + 7;
        auto rr = [&]() { s = s * 1103515245u + 12345u; return (int)((s >> 16) & 0x7fff); };
        for (int t = 0; t < 200; t++) {
            int i = rr() % GW, j = rr() % GH; if (!inMask[idx(i, j)]) continue;
            bool edge = false;
            for (int dj = -1; dj <= 1 && !edge; dj++) for (int di = -1; di <= 1; di++) {
                int ni = i + di, nj = j + dj; if (ni < 0 || ni >= GW || nj < 0 || nj >= GH || !inMask[idx(ni, nj)]) { edge = true; break; }
            }
            if (edge) { best = idx(i, j); break; }
        }
        if (best >= 0) heat[best] = 1.4f;
    }

    void step(real dt, const Vector3& wind) {
        time += dt;
        if (!ignited && (time >= igniteAt || pos.y < igniteY)) { igniteEdge(); ignited = true; }
        aero(dt, wind);
        if (ignited) fire(dt);
        reframe(dt);
    }

    void aero(real dt, const Vector3& wind) {
        Vector3 g(0, -9.81, 0);
        Vector3 va = vel - wind; real speed = va.magnitude();
        Vector3 n = e2; real vn = va * n; Vector3 vt = va - n * vn;
        real area = leafArea * (0.5 + 0.5 * aliveFrac());          // burnt leaves keep some drag
        Vector3 F = g * mass;
        F -= n * (kNormal * area * std::fabs(vn) * vn);               // face drag
        F -= vt * (kTang * area * vt.magnitude());                   // edge drag
        if (speed > 1e-4) {                                          // reactive lift ⟂ airflow
            Vector3 vah = va * (1.0 / speed);
            Vector3 lift = (n % vah) % vah; real lm = lift.magnitude();
            if (lm > 1e-6) F += lift * (kLift * area * speed * vn / lm);
        }
        vel += F * (dt / mass); vel *= airDamp;
        pos += vel * dt;
        if (speed > 1e-4) {                                          // torque turns the leaf broadside, underdamped
            Vector3 vah = va * (1.0 / speed);
            Vector3 tau = (n % (-vah)) * (kAlign * area * speed * speed) - omega * kAngDamp;
            omega += tau * (dt / inertia);
        }
    }

    void reframe(real dt) {
        e0 += (omega % e0) * dt; e1 += (omega % e1) * dt;
        e0.normalise(); e1 -= e0 * (e1 * e0); e1.normalise(); e2 = (e0 % e1);
    }

    void fire(real dt) {
        std::vector<float> add(GW * GH, 0.0f);
        aliveCount = 0;
        for (int j = 0; j < GH; j++) for (int i = 0; i < GW; i++) {
            int c = idx(i, j); if (!inMask[c]) continue; if (fuel[c] > 0.04f) aliveCount++;
            if (fuel[c] <= 0) continue;
            if (heat[c] > 0.12f) charAmt[c] = std::min(1.0f, charAmt[c] + (float)dt * heat[c] * 3.5f);
            curl[c] = charAmt[c] * (float)curlAmt * (float)ch * GH * 0.5f;
            if (heat[c] > ignite) {
                fuel[c] -= (float)dt * (float)burnRate * (0.5f + heat[c]) / resist[c];
                if (fuel[c] <= 0) { fuel[c] = 0; continue; }
                for (int dj = -1; dj <= 1; dj++) for (int di = -1; di <= 1; di++) {
                    if (!di && !dj) continue; int ni = i + di, nj = j + dj;
                    if (ni < 0 || ni >= GW || nj < 0 || nj >= GH) continue; int k = idx(ni, nj);
                    if (!inMask[k] || fuel[k] <= 0) continue;
                    add[k] += (float)dt * (float)spread * heat[c] / (di && dj ? 1.41f : 1.0f) / resist[k];   // ragged, undulating front
                }
            }
        }
        for (int c = 0; c < GW * GH; c++) heat[c] = (heat[c] + add[c]) * (1.0f - (float)dt * (float)cool);
    }

    // world position of grid node (i,j), including out-of-plane curl
    Vector3 node(int i, int j) const {
        real lx = (i - (GW - 1) * 0.5) * cw, ly = (j - (GH - 1) * 0.5) * ch, lz = curl[idx(i, j)];
        return pos + e0 * lx + e1 * ly + e2 * lz;
    }
    bool isBurning(int c) const { return inMask[c] && fuel[c] > 0.04f && fuel[c] < 0.85f && heat[c] > ignite * 0.8f; }
    bool consumed() const { return aliveFrac() < 0.06; }
};

} // namespace phys
