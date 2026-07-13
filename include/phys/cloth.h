// Position-based cloth (Verlet integration + iterative distance constraints).
// Built on the engine's Vector3. Structural, shear and bend constraints give it
// realistic drape; it collides with spheres and a ground plane and responds to
// wind. Per-vertex normals are recomputed each step for shading.
#pragma once
#include "core.h"
#include <vector>

namespace phys {

class Cloth {
public:
    int W, H;                                  // particle grid dimensions
    std::vector<Vector3> pos, prev, normal;
    std::vector<unsigned char> pinned;
    struct Constraint { int a, b; real rest; };
    std::vector<Constraint> constraints;
    std::vector<unsigned> tris;                // triangle indices for rendering

    Vector3 gravity = Vector3(0, -9.81, 0);
    real damping = (real)0.992;
    int iterations = 24;

    struct Sphere { Vector3 c; real r; };
    std::vector<Sphere> spheres;
    real groundY = -1e18;
    real collisionMargin = (real)0.02;
    real friction = (real)0.45;                // tangential grip at contacts

    int idx(int x, int y) const { return y * W + x; }

    // Build a W×H grid in the plane spanned by (right, down) from `origin`.
    Cloth(int W, int H, real spacing, const Vector3& origin,
          const Vector3& right = Vector3(1, 0, 0), const Vector3& down = Vector3(0, -1, 0))
        : W(W), H(H) {
        pos.resize(W * H); prev.resize(W * H); normal.resize(W * H); pinned.assign(W * H, 0);
        for (int y = 0; y < H; y++) for (int x = 0; x < W; x++) {
            Vector3 p = origin + right * (spacing * x) + down * (spacing * y);
            pos[idx(x, y)] = prev[idx(x, y)] = p;
        }
        auto add = [&](int a, int b) { constraints.push_back({a, b, (pos[a] - pos[b]).magnitude()}); };
        for (int y = 0; y < H; y++) for (int x = 0; x < W; x++) {
            if (x + 1 < W) add(idx(x, y), idx(x + 1, y));            // structural
            if (y + 1 < H) add(idx(x, y), idx(x, y + 1));
            if (x + 1 < W && y + 1 < H) { add(idx(x, y), idx(x + 1, y + 1)); add(idx(x + 1, y), idx(x, y + 1)); } // shear
            if (x + 2 < W) add(idx(x, y), idx(x + 2, y));            // bend
            if (y + 2 < H) add(idx(x, y), idx(x, y + 2));
        }
        for (int y = 0; y + 1 < H; y++) for (int x = 0; x + 1 < W; x++) {
            unsigned a = idx(x, y), b = idx(x + 1, y), c = idx(x, y + 1), d = idx(x + 1, y + 1);
            tris.insert(tris.end(), {a, c, b, b, c, d});
        }
        computeNormals();
    }

    void pin(int x, int y) { pinned[idx(x, y)] = 1; }
    void unpin(int x, int y) { pinned[idx(x, y)] = 0; }
    void setPin(int x, int y, const Vector3& p) { int i = idx(x, y); pos[i] = prev[i] = p; pinned[i] = 1; }

    void step(real dt, const Vector3& wind = Vector3()) {
        // Verlet integrate with wind biased by how much each vertex faces it
        for (int i = 0; i < (int)pos.size(); i++) {
            if (pinned[i]) { prev[i] = pos[i]; continue; }
            Vector3 acc = gravity;
            real face = real_abs(normal[i] * wind.unit());
            acc += wind * ((real)0.4 + (real)0.6 * face);
            Vector3 temp = pos[i];
            pos[i] += (pos[i] - prev[i]) * damping + acc * (dt * dt);
            prev[i] = temp;
        }
        for (int it = 0; it < iterations; it++) {
            for (auto& c : constraints) {
                Vector3 delta = pos[c.b] - pos[c.a];
                real len = delta.magnitude(); if (len < 1e-9) continue;
                real diff = (len - c.rest) / len;
                bool pa = pinned[c.a], pb = pinned[c.b];
                if (pa && pb) continue;
                if (!pa && !pb) { pos[c.a] += delta * (diff * (real)0.5); pos[c.b] -= delta * (diff * (real)0.5); }
                else if (pa) pos[c.b] -= delta * diff;
                else pos[c.a] += delta * diff;
            }
            collide();
        }
        applyContactFriction();
        computeNormals();
    }
    // Damp tangential velocity where the cloth touches an obstacle, so it grips
    // rather than sliding forever off a frictionless surface.
    void applyContactFriction() {
        for (int i = 0; i < (int)pos.size(); i++) {
            if (pinned[i]) continue;
            Vector3 nrm; bool contact = false;
            for (auto& s : spheres) { Vector3 d = pos[i] - s.c; real len = d.magnitude();
                if (len < s.r + collisionMargin * 2 && len > 1e-9) { nrm = d * (((real)1) / len); contact = true; } }
            if (pos[i].y < groundY + collisionMargin * 2) { nrm = Vector3(0, 1, 0); contact = true; }
            if (!contact) continue;
            Vector3 vel = pos[i] - prev[i];
            Vector3 vn = nrm * (vel * nrm);
            Vector3 vt = vel - vn;
            prev[i] = pos[i] - (vn + vt * (1 - friction));
        }
    }
    void collide() {
        for (int i = 0; i < (int)pos.size(); i++) {
            if (pinned[i]) continue;
            for (auto& s : spheres) {
                Vector3 d = pos[i] - s.c; real len = d.magnitude(); real r = s.r + collisionMargin;
                if (len < r && len > 1e-9) pos[i] = s.c + d * (r / len);
            }
            if (pos[i].y < groundY + collisionMargin) pos[i].y = groundY + collisionMargin;
        }
    }
    void computeNormals() {
        for (auto& n : normal) n = Vector3();
        for (size_t t = 0; t < tris.size(); t += 3) {
            int a = tris[t], b = tris[t + 1], c = tris[t + 2];
            Vector3 fn = (pos[b] - pos[a]) % (pos[c] - pos[a]);
            normal[a] += fn; normal[b] += fn; normal[c] += fn;
        }
        for (auto& n : normal) { if (n.squareMagnitude() > 1e-12) n.normalise(); else n = Vector3(0, 0, 1); }
    }
};

} // namespace phys
