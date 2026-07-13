// Burning paper: a Verlet cloth sheet coupled to a fire-propagation cellular
// automaton. Each cell has fuel, heat and accumulated char. Heat diffuses (with
// an upward bias — fire climbs), ignites neighbouring cells, and consumes fuel;
// a consumed cell breaks its constraints (a hole opens and charred scraps fall).
// Charring cells curl out of plane and hot cells get a little buoyant lift.
#pragma once
#include "core.h"
#include <vector>
#include <cmath>

namespace phys {

class BurningPaper {
public:
    int W, H;
    std::vector<Vector3> pos, prev, normal;
    std::vector<unsigned char> pinned;
    std::vector<float> fuel;        // 1 = fresh paper, 0 = consumed (hole)
    std::vector<float> heat;        // temperature field
    std::vector<float> charAmt;     // 0 = white, 1 = black char (permanent)
    std::vector<float> resist;      // per-cell burn resistance → ragged front
    struct Constraint { int a, b; float rest; bool active; };
    std::vector<Constraint> constraints;
    std::vector<unsigned> tris;     // quad triangles (drawn only where alive)

    Vector3 gravity = Vector3(0, -3.2, 0);   // paper is light
    int iterations = 14;
    float damping = 0.985f;
    // fire tuning
    float ignite = 0.5f, burnRate = 0.62f, spread = 3.4f, upBias = 2.4f;
    float cool = 0.55f, buoyancy = 4.5f, curlForce = 9.0f;
    float flutter = 0.0f, flutterSpeed = 3.5f;   // gentle out-of-plane breeze that ripples across the sheet
    double time = 0;

    int idx(int x, int y) const { return y * W + x; }
    bool cellAlive(int i) const { return fuel[i] > 0.04f; }

    void build(int W_, int H_, real spacing, const Vector3& origin,
               const Vector3& right = Vector3(1, 0, 0), const Vector3& down = Vector3(0, -1, 0)) {
        W = W_; H = H_;
        int n = W * H;
        pos.resize(n); prev.resize(n); normal.resize(n); pinned.assign(n, 0);
        fuel.assign(n, 1); heat.assign(n, 0); charAmt.assign(n, 0); resist.resize(n);
        unsigned s = 1234567u; auto rr = [&]() { s = s * 1103515245u + 12345u; return (float)(((s >> 16) & 0x7fff) / 32767.0); };
        for (int y = 0; y < H; y++) for (int x = 0; x < W; x++) {
            pos[idx(x, y)] = prev[idx(x, y)] = origin + right * (spacing * x) + down * (spacing * y);
            resist[idx(x, y)] = 0.7f + 0.6f * rr();
        }
        auto add = [&](int a, int b) { constraints.push_back({a, b, (float)(pos[a] - pos[b]).magnitude(), true}); };
        for (int y = 0; y < H; y++) for (int x = 0; x < W; x++) {
            if (x + 1 < W) add(idx(x, y), idx(x + 1, y));
            if (y + 1 < H) add(idx(x, y), idx(x, y + 1));
            if (x + 1 < W && y + 1 < H) { add(idx(x, y), idx(x + 1, y + 1)); add(idx(x + 1, y), idx(x, y + 1)); }
            if (x + 2 < W) add(idx(x, y), idx(x + 2, y));
            if (y + 2 < H) add(idx(x, y), idx(x, y + 2));
        }
        for (int y = 0; y + 1 < H; y++) for (int x = 0; x + 1 < W; x++) {
            unsigned a = idx(x, y), b = idx(x + 1, y), c = idx(x, y + 1), d = idx(x + 1, y + 1);
            tris.insert(tris.end(), {a, c, b, b, c, d});
        }
        computeNormals();
    }
    void pin(int x, int y) { pinned[idx(x, y)] = 1; }
    void ignite_at(int x, int y, float amount) { if (x >= 0 && x < W && y >= 0 && y < H) heat[idx(x, y)] += amount; }

    void step(real dt, const Vector3& wind = Vector3()) {
        time += dt;
        // ---- fire cellular automaton ----
        std::vector<float> add(W * H, 0.0f);
        for (int y = 0; y < H; y++) for (int x = 0; x < W; x++) {
            int i = idx(x, y); if (fuel[i] <= 0) continue;
            bool burning = heat[i] > ignite && fuel[i] > 0;
            if (heat[i] > 0.14f) charAmt[i] = std::min(1.0f, charAmt[i] + (float)dt * heat[i] * 1.8f);
            if (burning) {
                fuel[i] -= (float)dt * burnRate * (0.5f + 0.9f * heat[i]) / resist[i];
                if (fuel[i] <= 0) { fuel[i] = 0; deactivateAround(i); continue; }
                // emit heat to 8 neighbours, biased upward (row-1 = up in world)
                for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
                    if (!dx && !dy) continue; int nx = x + dx, ny = y + dy;
                    if (nx < 0 || nx >= W || ny < 0 || ny >= H) continue; int j = idx(nx, ny);
                    if (fuel[j] <= 0) continue;
                    float w = (dy < 0 ? upBias : (dy > 0 ? 0.5f : 1.0f)) / (dx && dy ? 1.41f : 1.0f);
                    add[j] += (float)dt * spread * heat[i] * w / resist[j];
                }
            }
        }
        for (int i = 0; i < W * H; i++) heat[i] = (heat[i] + add[i]) * (1.0f - (float)dt * cool);

        // ---- cloth physics ----
        for (int i = 0; i < W * H; i++) {
            if (pinned[i] || fuel[i] <= 0) { prev[i] = pos[i]; continue; }
            Vector3 acc = gravity + wind;
            if (heat[i] > 0.1f) acc.y += buoyancy * heat[i];                    // hot air lifts
            if (flutter > 0.0f) {                                              // travelling ripple, stronger toward the free lower edge
                int x = i % W, y = i / W;
                float phase = (float)x * 0.42f + (float)y * 0.26f - (float)time * flutterSpeed;
                acc.z += flutter * std::sin(phase) * (0.35f + (float)y / (float)H);
            }
            bool front = fuel[i] < 0.75f && heat[i] > ignite * 0.7f;
            if (front) acc += normal[i] * (curlForce * (0.75f - fuel[i]));      // char curls out of plane
            Vector3 tmp = pos[i];
            pos[i] += (pos[i] - prev[i]) * damping + acc * (real)(dt * dt);
            prev[i] = tmp;
        }
        for (int it = 0; it < iterations; it++)
            for (auto& c : constraints) {
                if (!c.active) continue;
                Vector3 d = pos[c.b] - pos[c.a]; real len = d.magnitude(); if (len < 1e-9) continue;
                real diff = (len - c.rest) / len;
                bool pa = pinned[c.a], pb = pinned[c.b];
                if (pa && pb) continue;
                if (!pa && !pb) { pos[c.a] += d * (diff * (real)0.5); pos[c.b] -= d * (diff * (real)0.5); }
                else if (pa) pos[c.b] -= d * diff;
                else pos[c.a] += d * diff;
            }
        computeNormals();
    }
    void deactivateAround(int cell) {
        for (auto& c : constraints) if (c.a == cell || c.b == cell) c.active = false;
    }
    void computeNormals() {
        for (auto& n : normal) n = Vector3();
        for (size_t t = 0; t + 2 < tris.size(); t += 3) {
            int a = tris[t], b = tris[t + 1], c = tris[t + 2];
            if (fuel[a] <= 0.04f || fuel[b] <= 0.04f || fuel[c] <= 0.04f) continue;
            Vector3 fn = (pos[b] - pos[a]) % (pos[c] - pos[a]);
            normal[a] += fn; normal[b] += fn; normal[c] += fn;
        }
        for (auto& n : normal) { if (n.squareMagnitude() > 1e-12) n.normalise(); else n = Vector3(0, 0, 1); }
    }
    // for spawning smoke/embers: is this cell actively on fire?
    bool isBurning(int i) const { return fuel[i] > 0.04f && fuel[i] < 0.9f && heat[i] > ignite * 0.8f; }
};

} // namespace phys
