// Rope simulation: three ropes doing rope things — (1) a rope pinned to a beam
// swinging a heavy ball, (2) a rope draped over the beam with both ends hanging,
// (3) a rope dropped from height that coils on the floor. Verlet chains with
// inextensible segments, ground + beam collision and tangential friction.
//   ./rope3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
using namespace phys;

struct Rope {
    std::vector<Vector3> pos, prev;
    std::vector<unsigned char> pinned;
    real seg, radius; real tipMass = 1;             // mass ratio of the last particle
    void build(const Vector3& a, const Vector3& b, int n, real r) {
        pos.clear(); prev.clear(); pinned.clear(); radius = r;
        for (int i = 0; i <= n; i++) { Vector3 p = a + (b - a) * ((real)i / n); pos.push_back(p); prev.push_back(p); pinned.push_back(0); }
        seg = (b - a).magnitude() / n;
    }
    void step(real dt, const Vector3& beamC, real beamR) {
        for (int i = 0; i < (int)pos.size(); i++) {
            if (pinned[i]) { prev[i] = pos[i]; continue; }
            Vector3 tmp = pos[i];
            pos[i] += (pos[i] - prev[i]) * (real)0.985 + Vector3(0, -9.81, 0) * (dt * dt);
            prev[i] = tmp;
        }
        for (int it = 0; it < 30; it++) {
            for (int i = 0; i + 1 < (int)pos.size(); i++) {           // inextensible segments
                Vector3 d = pos[i + 1] - pos[i]; real len = d.magnitude(); if (len < 1e-9) continue;
                real w0 = pinned[i] ? 0 : 1, w1 = pinned[i + 1] ? 0 : (i + 2 == (int)pos.size() ? 1 / tipMass : 1);
                real tw = w0 + w1; if (tw <= 0) continue;
                Vector3 corr = d * ((len - seg) / len);
                pos[i] += corr * (w0 / tw); pos[i + 1] -= corr * (w1 / tw);
            }
            for (int i = 0; i < (int)pos.size(); i++) {               // collisions
                if (pinned[i]) continue;
                if (pos[i].y < radius) {                              // ground + friction
                    pos[i].y = radius;
                    Vector3 v = pos[i] - prev[i]; v.y = 0; prev[i] = pos[i] - v * (real)0.4;
                }
                Vector3 d = pos[i] - beamC; d.x = 0;                  // beam runs along x
                real len = d.magnitude();
                if (len < beamR + radius && len > 1e-9) pos[i] = Vector3(pos[i].x, beamC.y, beamC.z) + d * ((beamR + radius) / len);
            }
        }
    }
};

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const Vector3 beamC(0, 3.0, 0); const real beamR = 0.09;

    Rope pendulum;                                    // swings a heavy ball
    pendulum.build(Vector3(-1.6, 3.0, 0), Vector3(-3.4, 2.2, 1.6), 26, 0.042);
    pendulum.pinned[0] = 1; pendulum.tipMass = 10;
    Rope draped;                                      // hangs over the beam
    draped.build(Vector3(0.8, 3.35, -1.5), Vector3(0.8, 3.35, 1.5), 34, 0.042);
    Rope faller;                                      // coils on the floor
    faller.build(Vector3(2.6, 4.6, -0.4), Vector3(2.6, 1.2, 0.5), 34, 0.042);

    gfx::App app(W, H, "ropes", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(30, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.7, 0); cam.dist = 8.2f; cam.yaw = 0.5f; cam.pitch = 0.2f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    // --- realistic rope: swept tube with a twisted 3-strand cross-section ---
    struct RopeMesh {
        gfx::Mesh mesh; std::vector<float> verts; std::vector<unsigned> idx; bool init = false;
        static Vector3 catmull(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, real t) {
            real t2 = t * t, t3 = t2 * t;
            return (p1 * 2 + (p2 - p0) * t + (p0 * 2 - p1 * 5 + p2 * 4 - p3) * t2 + (p1 * 3 - p0 - p2 * 3 + p3) * t3) * (real)0.5;
        }
        void build(const Rope& rp) {
            const int SUB = 3, RING = 14; const real TWIST = 34.0;   // strand turns per unit length
            std::vector<Vector3> c;                                    // smoothed centreline
            int n = (int)rp.pos.size();
            for (int i = 0; i + 1 < n; i++) {
                const Vector3& p0 = rp.pos[i > 0 ? i - 1 : 0]; const Vector3& p1 = rp.pos[i];
                const Vector3& p2 = rp.pos[i + 1]; const Vector3& p3 = rp.pos[i + 2 < n ? i + 2 : n - 1];
                for (int s = 0; s < SUB; s++) c.push_back(catmull(p0, p1, p2, p3, (real)s / SUB));
            }
            c.push_back(rp.pos.back());
            int R = (int)c.size();
            // parallel-transported frame
            Vector3 t0 = (c[1] - c[0]).unit();
            Vector3 up = real_abs(t0.y) < 0.95 ? Vector3(0, 1, 0) : Vector3(1, 0, 0);
            Vector3 nrm = (t0 % up).unit();
            verts.assign((size_t)R * RING * 6, 0.0f);
            real arc = 0;
            std::vector<Vector3> grid((size_t)R * RING);
            for (int i = 0; i < R; i++) {
                Vector3 tan = (c[i < R - 1 ? i + 1 : i] - c[i > 0 ? i - 1 : i]).unit();
                nrm = (nrm - tan * (nrm * tan)).unit();               // transport
                Vector3 bin = (tan % nrm).unit();
                if (i > 0) arc += (c[i] - c[i - 1]).magnitude();
                for (int k = 0; k < RING; k++) {
                    real th = 2 * real_pi * k / RING;
                    real bump = 1 + (real)0.16 * real_cos(3 * th - TWIST * arc);   // 3-strand helix
                    real rr2 = rp.radius * bump;
                    grid[(size_t)i * RING + k] = c[i] + (nrm * real_cos(th) + bin * real_sin(th)) * rr2;
                }
            }
            idx.clear();
            for (int i = 0; i + 1 < R; i++) for (int k = 0; k < RING; k++) {
                unsigned a = i * RING + k, b = i * RING + (k + 1) % RING, d = (i + 1) * RING + k, e = (i + 1) * RING + (k + 1) % RING;
                idx.insert(idx.end(), {a, d, b, b, d, e});
            }
            std::vector<Vector3> vn((size_t)R * RING);
            for (size_t f = 0; f < idx.size(); f += 3) {
                Vector3 fn = (grid[idx[f + 1]] - grid[idx[f]]) % (grid[idx[f + 2]] - grid[idx[f]]);
                vn[idx[f]] += fn; vn[idx[f + 1]] += fn; vn[idx[f + 2]] += fn;
            }
            for (size_t v = 0; v < grid.size(); v++) { Vector3 nn = vn[v]; if (nn.squareMagnitude() > 1e-12) nn.normalise();
                float* o = &verts[v * 6];
                o[0] = grid[v].x; o[1] = grid[v].y; o[2] = grid[v].z; o[3] = nn.x; o[4] = nn.y; o[5] = nn.z; }
            if (!init) { mesh.upload(verts, idx); init = true; } else mesh.updateVerts(verts);
        }
    };
    RopeMesh rmPend, rmDrape, rmFall;
    auto drawRope = [&](RopeMesh& rm, glm::vec3 col, bool shadow) {
        if (shadow) r.shadowDraw(rm.mesh, glm::mat4(1));
        else r.drawPBR(r.pSolid, rm.mesh, glm::mat4(1), col, 0.92f, 0.0f);
    };
    glm::mat4 beamM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(beamC.x, beamC.y, beamC.z)), glm::vec3(4.0f, beamR, beamR));
    glm::mat4 postL = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-4.0f, 1.5f, 0)), glm::vec3(0.1f, 1.5f, 0.1f));
    glm::mat4 postR = glm::scale(glm::translate(glm::mat4(1), glm::vec3(4.0f, 1.5f, 0)), glm::vec3(0.1f, 1.5f, 0.1f));

    auto renderFrame = [&]() {
        Vector3 ballP = pendulum.pos.back();
        glm::mat4 ballM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(ballP.x, ballP.y, ballP.z)), glm::vec3(0.22f));
        r.setLightForScene(glm::vec3(0, 1.6, 0), 7.0f);
        rmPend.build(pendulum); rmDrape.build(draped); rmFall.build(faller);
        r.beginShadow();
        r.shadowDraw(box, beamM); r.shadowDraw(box, postL); r.shadowDraw(box, postR); r.shadowDraw(sphere, ballM);
        drawRope(rmPend, {}, true); drawRope(rmDrape, {}, true); drawRope(rmFall, {}, true);
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.58, 0.56, 0.52), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, beamM, glm::vec3(0.4, 0.28, 0.16), 0.8f, 0.0f);
        r.drawPBR(r.pSolid, box, postL, glm::vec3(0.4, 0.28, 0.16), 0.8f, 0.0f);
        r.drawPBR(r.pSolid, box, postR, glm::vec3(0.4, 0.28, 0.16), 0.8f, 0.0f);
        r.drawPBR(r.pSolid, sphere, ballM, glm::vec3(0.25, 0.27, 0.32), 0.3f, 0.8f);
        drawRope(rmPend, glm::vec3(0.72, 0.58, 0.34), false);          // hemp
        drawRope(rmDrape, glm::vec3(0.66, 0.52, 0.30), false);
        drawRope(rmFall, glm::vec3(0.58, 0.46, 0.27), false);
        r.endScene();
    };
    auto stepAll = [&](real dt) { pendulum.step(dt, beamC, beamR); draped.step(dt, beamC, beamR); faller.step(dt, beamC, beamR); };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 0.5f - f * (0.3f / frames); renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 4; s++) stepAll(1.0 / 240); }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames * 4; i++) stepAll(1.0 / 240); renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        for (int s = 0; s < 4; s++) stepAll(1.0 / 240);
        renderFrame(); r.present(); app.poll();
    }
    return 0;
}
