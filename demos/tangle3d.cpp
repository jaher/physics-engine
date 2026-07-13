// Chains hanging and getting entangled: four heavy iron chains hang from anchors
// that orbit a common centre (a maypole braid). Full chain-vs-chain collision —
// every link segment repels every other with thickness + friction — makes them
// catch, wrap and knot around each other, then hold the tangle when the motion
// stops. Links drawn as interlocked oval tori.
//   ./tangle3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
#include <vector>
using namespace phys;

// closest points between two segments (clamped) — for chain-chain collision
static void closestSeg(const Vector3& p1, const Vector3& q1, const Vector3& p2, const Vector3& q2,
                       Vector3& c1, Vector3& c2) {
    Vector3 d1 = q1 - p1, d2 = q2 - p2, rr = p1 - p2;
    real a = d1 * d1, e = d2 * d2, f = d2 * rr, s, t;
    if (a <= 1e-9 && e <= 1e-9) { c1 = p1; c2 = p2; return; }
    if (a <= 1e-9) { s = 0; t = f / e; }
    else { real c = d1 * rr;
        if (e <= 1e-9) { t = 0; s = -c / a; }
        else { real b = d1 * d2, den = a * e - b * b;
            s = den > 1e-12 ? (b * f - c * e) / den : 0;
            s = s < 0 ? 0 : (s > 1 ? 1 : s);
            t = (b * s + f) / e;
            if (t < 0) { t = 0; s = -c / a; } else if (t > 1) { t = 1; s = (b - c) / a; } } }
    s = s < 0 ? 0 : (s > 1 ? 1 : s); c1 = p1 + d1 * s; c2 = p2 + d2 * t;
}

struct Chain {
    std::vector<Vector3> pos, prev; std::vector<unsigned char> pinned;
    real pitch, thick;
    void build(const Vector3& a, const Vector3& b, int n, real th) {
        thick = th; pos.clear(); prev.clear(); pinned.clear();
        for (int i = 0; i <= n; i++) { Vector3 p = a + (b - a) * ((real)i / n); pos.push_back(p); prev.push_back(p); pinned.push_back(0); }
        pitch = (b - a).magnitude() / n;
    }
    void integrate(real dt) {
        for (int i = 0; i < (int)pos.size(); i++) {
            if (pinned[i]) { prev[i] = pos[i]; continue; }
            Vector3 tmp = pos[i];
            pos[i] += (pos[i] - prev[i]) * (real)0.99 + Vector3(0, -9.81, 0) * (dt * dt);
            prev[i] = tmp;
        }
    }
    void satisfy() {
        for (int i = 0; i + 1 < (int)pos.size(); i++) {
            Vector3 d = pos[i + 1] - pos[i]; real len = d.magnitude(); if (len < 1e-9) continue;
            real w0 = pinned[i] ? 0 : 1, w1 = pinned[i + 1] ? 0 : 1; real tw = w0 + w1; if (tw <= 0) continue;
            Vector3 corr = d * ((len - pitch) / len);
            pos[i] += corr * (w0 / tw); pos[i + 1] -= corr * (w1 / tw);
        }
    }
};

// push two capsule-segments apart to their combined thickness, with friction
static void collideSeg(Chain& A, int i, Chain& B, int j, real dt) {
    Vector3 a0 = A.pos[i], a1 = A.pos[i + 1], b0 = B.pos[j], b1 = B.pos[j + 1];
    Vector3 ca, cb; closestSeg(a0, a1, b0, b1, ca, cb);
    Vector3 d = ca - cb; real len = d.magnitude();
    real minDist = A.thick + B.thick;
    if (len >= minDist || len < 1e-9) return;
    Vector3 n = d * (((real)1) / len);
    real pen = minDist - len;
    // distribute the correction to the four endpoints by proximity to the contact
    real ta = (ca - a0).magnitude() / (A.pitch + 1e-9); ta = ta < 0 ? 0 : (ta > 1 ? 1 : ta);
    real tb = (cb - b0).magnitude() / (B.pitch + 1e-9); tb = tb < 0 ? 0 : (tb > 1 ? 1 : tb);
    Vector3 push = n * (pen * (real)0.22);   // soft push → smooth wrap, no jitter
    auto add = [](Chain& C, int idx, const Vector3& v) { if (!C.pinned[idx]) C.pos[idx] += v; };
    add(A, i, push * (1 - ta)); add(A, i + 1, push * ta);
    add(B, j, push * -(1 - tb)); add(B, j + 1, push * -tb);
    // tangential friction: damp the sliding of the contact points
    Vector3 va = ((A.pos[i] - A.prev[i]) * (1 - ta) + (A.pos[i + 1] - A.prev[i + 1]) * ta);
    Vector3 vb = ((B.pos[j] - B.prev[j]) * (1 - tb) + (B.pos[j + 1] - B.prev[j + 1]) * tb);
    Vector3 rel = va - vb; Vector3 tang = rel - n * (rel * n);
    Vector3 fr = tang * (real)0.12;
    if (!A.pinned[i]) A.prev[i] += fr * (1 - ta); if (!A.pinned[i + 1]) A.prev[i + 1] += fr * ta;
    if (!B.pinned[j]) B.prev[j] -= fr * (1 - tb); if (!B.pinned[j + 1]) B.prev[j + 1] -= fr * tb;
    (void)dt;
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 360;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int NLINK = 44; const real THICK = 0.05, TOP = 3.4, RING = 0.9, YBOT = 0.7, SLACK = 1.14;

    std::vector<Chain> chains(4);
    glm::vec3 cols[4] = {{0.33,0.32,0.34},{0.35,0.29,0.25},{0.27,0.30,0.35},{0.34,0.33,0.28}};
    for (int c = 0; c < 4; c++) {
        double a = c * glm::half_pi<double>();
        Vector3 top(std::cos(a) * RING, TOP, std::sin(a) * RING);
        Vector3 bot(std::cos(a) * RING, YBOT, std::sin(a) * RING);
        chains[c].build(top, bot, NLINK, THICK);
        chains[c].pitch *= SLACK;                         // slack → gentle belly, room to wrap
        chains[c].pinned[0] = 1; chains[c].pinned[NLINK] = 1;
    }

    gfx::App app(W, H, "entangling chains", headless);
    gfx::Renderer r; r.init(W, H);
    r.lightDir = glm::normalize(glm::vec3(-0.6f, -0.8f, -0.4f));
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(30, 1);
    gfx::Mesh torus = gfx::makeTorus(1.0f, 0.30f, 32, 16);
    glm::mat4 hub = glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, TOP + 0.12f, 0)), glm::vec3(1.15f, 0.1f, 1.15f));

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 2.0, 0); cam.dist = 6.8f; cam.yaw = 0.5f; cam.pitch = 0.14f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    double t = 0; double curTwist = 0;
    auto stepPhysics = [&](real dt) {
        double twist = (t < 15.0 ? t : 15.0) * 0.62;      // bottom ring winds ~1.5 turns, then holds
        curTwist = twist;
        for (int c = 0; c < 4; c++) {
            double a = c * glm::half_pi<double>();
            chains[c].pos[0] = Vector3(std::cos(a) * RING, TOP, std::sin(a) * RING);
            chains[c].pos[NLINK] = Vector3(std::cos(a + twist) * RING, YBOT, std::sin(a + twist) * RING);
            chains[c].integrate(dt);
        }
        for (int it = 0; it < 22; it++) {
            for (auto& ch : chains) ch.satisfy();
            for (auto& ch : chains) for (int i = 0; i < (int)ch.pos.size(); i++)
                if (!ch.pinned[i] && ch.pos[i].y < THICK) { ch.pos[i].y = THICK;
                    Vector3 v = ch.pos[i] - ch.prev[i]; v.y = 0; ch.prev[i] = ch.pos[i] - v * (real)0.5; }
            // chain-chain collision (all pairs, skip adjacent segments on the same chain)
            for (int c1 = 0; c1 < 4; c1++) for (int c2 = c1; c2 < 4; c2++)
                for (int i = 0; i + 1 < (int)chains[c1].pos.size(); i++)
                    for (int j = (c1 == c2 ? i + 2 : 0); j + 1 < (int)chains[c2].pos.size(); j++)
                        collideSeg(chains[c1], i, chains[c2], j, dt);
        }
        t += dt;
    };

    auto linkMat = [&](const Chain& ch, int i) {
        Vector3 a = ch.pos[i], b = ch.pos[i + 1];
        Vector3 mid = (a + b) * 0.5, d = b - a; real len = d.magnitude(); if (len < 1e-9) d = Vector3(0, 1, 0); else d.normalise();
        glm::vec3 dir((float)d.x, (float)d.y, (float)d.z);
        glm::vec3 up = std::fabs(dir.y) < 0.95f ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
        glm::vec3 x = glm::normalize(glm::cross(up, dir)), y = glm::cross(dir, x);
        glm::mat4 M(glm::vec4(x, 0), glm::vec4(y, 0), glm::vec4(dir, 0), glm::vec4(mid.x, mid.y, mid.z, 1));
        M = glm::rotate(M, (float)(i % 2) * glm::half_pi<float>(), glm::vec3(0, 0, 1));
        M = glm::rotate(M, glm::half_pi<float>(), glm::vec3(1, 0, 0));
        return glm::scale(M, glm::vec3((float)ch.pitch * 0.5f, (float)ch.pitch * 0.72f, (float)ch.pitch * 0.5f));
    };
    auto drawChains = [&](bool shadow) {
        for (int c = 0; c < 4; c++) for (int i = 0; i + 1 < (int)chains[c].pos.size(); i++) {
            glm::mat4 m = linkMat(chains[c], i);
            if (shadow) r.shadowDraw(torus, m);
            else r.drawPBR(r.pSolid, torus, m, cols[c], 0.5f, 0.9f);
        }
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 1.6, 0), 5.5f);
        glm::mat4 botPlate = glm::scale(glm::rotate(glm::translate(glm::mat4(1), glm::vec3(0, YBOT - 0.1f, 0)), (float)curTwist, glm::vec3(0,1,0)), glm::vec3(1.15f, 0.1f, 1.15f));
        r.beginShadow(); r.shadowDraw(box, hub); r.shadowDraw(box, botPlate); drawChains(true); r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.5, 0.49, 0.46), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, hub, glm::vec3(0.26, 0.25, 0.27), 0.5f, 0.9f);
        r.drawPBR(r.pSolid, box, botPlate, glm::vec3(0.26, 0.25, 0.27), 0.5f, 0.9f);
        drawChains(false);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 0.35f + f * (0.6f / frames); renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 7; s++) stepPhysics(1.0 / 280); }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames * 5; i++) stepPhysics(1.0 / 300); renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        for (int s = 0; s < 5; s++) stepPhysics(1.0 / 300);
        renderFrame(); r.present(); app.poll();
    }
    return 0;
}
