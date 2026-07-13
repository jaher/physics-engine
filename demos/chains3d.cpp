// Medieval iron chains: interlocked oval links (elongated tori, alternating 90°
// roll along the run) on Verlet chain physics. Three set-pieces: a heavy chain
// slung between two stone pillars, a flail chain swinging an iron ball from a
// gibbet arm, and a chain spilling into a pile on the flagstones.
//   ./chains3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
using namespace phys;

struct Chain {
    std::vector<Vector3> pos, prev;
    std::vector<unsigned char> pinned;
    real pitch; real tipMass = 1;
    void build(const Vector3& a, const Vector3& b, int n) {
        pos.clear(); prev.clear(); pinned.clear();
        for (int i = 0; i <= n; i++) { Vector3 p = a + (b - a) * ((real)i / n); pos.push_back(p); prev.push_back(p); pinned.push_back(0); }
        pitch = (b - a).magnitude() / n;
    }
    void step(real dt) {
        for (int i = 0; i < (int)pos.size(); i++) {
            if (pinned[i]) { prev[i] = pos[i]; continue; }
            Vector3 tmp = pos[i];
            pos[i] += (pos[i] - prev[i]) * (real)0.985 + Vector3(0, -9.81, 0) * (dt * dt);
            prev[i] = tmp;
        }
        for (int it = 0; it < 30; it++) {
            for (int i = 0; i + 1 < (int)pos.size(); i++) {
                Vector3 d = pos[i + 1] - pos[i]; real len = d.magnitude(); if (len < 1e-9) continue;
                real w0 = pinned[i] ? 0 : 1, w1 = pinned[i + 1] ? 0 : (i + 2 == (int)pos.size() ? 1 / tipMass : 1);
                real tw = w0 + w1; if (tw <= 0) continue;
                Vector3 corr = d * ((len - pitch) / len);
                pos[i] += corr * (w0 / tw); pos[i + 1] -= corr * (w1 / tw);
            }
            for (int i = 0; i < (int)pos.size(); i++) {               // ground with iron-on-stone friction
                if (pinned[i]) continue;
                if (pos[i].y < 0.05) { pos[i].y = 0.05;
                    Vector3 v = pos[i] - prev[i]; v.y = 0; prev[i] = pos[i] - v * (real)0.35; }
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

    Chain slung;                                      // pillar to pillar, sagging
    slung.build(Vector3(-3.2, 2.6, -0.6), Vector3(0.9, 2.6, -0.6), 30);
    slung.pinned[0] = 1; slung.pinned[30] = 1;
    { for (int i = 1; i < 30; i++) slung.pos[i].y -= 0.001 * i * (30 - i); }   // pre-sag nudge
    Chain flail;                                      // gibbet arm + iron ball
    flail.build(Vector3(2.6, 3.4, 0.2), Vector3(1.4, 2.2, 1.4), 16);
    flail.pinned[0] = 1; flail.tipMass = 14;
    Chain spill;                                      // piling onto the floor
    spill.build(Vector3(-1.4, 4.8, 1.6), Vector3(-1.5, 0.8, 2.4), 34);

    gfx::App app(W, H, "medieval chains", headless);
    gfx::Renderer r; r.init(W, H);
    r.lightDir = glm::normalize(glm::vec3(-0.7f, -0.75f, -0.4f));
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(30, 1);
    gfx::Mesh torus = gfx::makeTorus(1.0f, 0.30f, 40, 20);

    gfx::OrbitCamera cam; cam.target = glm::vec3(-0.4, 1.6, 0.4); cam.dist = 7.6f; cam.yaw = 0.35f; cam.pitch = 0.2f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    // one interlocked link between consecutive chain nodes, alternating 90° roll
    auto linkMat = [&](const Chain& ch, int i) {
        Vector3 a = ch.pos[i], b = ch.pos[i + 1];
        Vector3 mid = (a + b) * 0.5, d = b - a; real len = d.magnitude(); d.normalise();
        glm::vec3 dir((float)d.x, (float)d.y, (float)d.z);
        glm::vec3 up = std::fabs(dir.y) < 0.95f ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
        glm::vec3 x = glm::normalize(glm::cross(up, dir)), y = glm::cross(dir, x);
        glm::mat4 M(glm::vec4(x, 0), glm::vec4(y, 0), glm::vec4(dir, 0), glm::vec4(mid.x, mid.y, mid.z, 1));
        M = glm::rotate(M, (float)(i % 2) * glm::half_pi<float>(), glm::vec3(0, 0, 1));   // interlock
        // oval link: torus lies in local XZ, elongated along the chain (Z)
        M = glm::rotate(M, glm::half_pi<float>(), glm::vec3(1, 0, 0));
        return glm::scale(M, glm::vec3((float)ch.pitch * 0.42f, (float)ch.pitch * 0.72f, (float)ch.pitch * 0.42f));
    };
    auto drawChain = [&](const Chain& ch, bool shadow) {
        for (int i = 0; i + 1 < (int)ch.pos.size(); i++) {
            glm::mat4 m = linkMat(ch, i);
            if (shadow) r.shadowDraw(torus, m);
            else { float j = 0.9f + 0.1f * (float)((i * 2654435761u >> 8 & 0xff) / 255.0);
                r.drawPBR(r.pSolid, torus, m, glm::vec3(0.30f, 0.29f, 0.30f) * j, 0.48f, 0.92f); }
        }
    };
    glm::mat4 pillarL = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-3.4f, 1.3f, -0.6f)), glm::vec3(0.34f, 1.3f, 0.34f));
    glm::mat4 pillarR = glm::scale(glm::translate(glm::mat4(1), glm::vec3(1.1f, 1.3f, -0.6f)), glm::vec3(0.34f, 1.3f, 0.34f));
    glm::mat4 capL = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-3.4f, 2.68f, -0.6f)), glm::vec3(0.45f, 0.09f, 0.45f));
    glm::mat4 capR = glm::scale(glm::translate(glm::mat4(1), glm::vec3(1.1f, 2.68f, -0.6f)), glm::vec3(0.45f, 0.09f, 0.45f));
    glm::mat4 gibbetPost = glm::scale(glm::translate(glm::mat4(1), glm::vec3(3.1f, 1.7f, 0.2f)), glm::vec3(0.09f, 1.7f, 0.09f));
    glm::mat4 gibbetArm = glm::scale(glm::translate(glm::mat4(1), glm::vec3(2.8f, 3.4f, 0.2f)), glm::vec3(0.4f, 0.07f, 0.07f));

    auto renderFrame = [&]() {
        Vector3 ballP = flail.pos.back();
        glm::mat4 ballM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(ballP.x, ballP.y - 0.1f, ballP.z)), glm::vec3(0.24f));
        r.setLightForScene(glm::vec3(0, 1.5, 0.5), 7.0f);
        r.beginShadow();
        r.shadowDraw(box, pillarL); r.shadowDraw(box, pillarR); r.shadowDraw(box, capL); r.shadowDraw(box, capR);
        r.shadowDraw(box, gibbetPost); r.shadowDraw(box, gibbetArm); r.shadowDraw(sphere, ballM);
        drawChain(slung, true); drawChain(flail, true); drawChain(spill, true);
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.52, 0.50, 0.47), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, pillarL, glm::vec3(0.55, 0.52, 0.47), 0.9f, 0.0f);
        r.drawPBR(r.pSolid, box, pillarR, glm::vec3(0.55, 0.52, 0.47), 0.9f, 0.0f);
        r.drawPBR(r.pSolid, box, capL, glm::vec3(0.48, 0.45, 0.41), 0.9f, 0.0f);
        r.drawPBR(r.pSolid, box, capR, glm::vec3(0.48, 0.45, 0.41), 0.9f, 0.0f);
        r.drawPBR(r.pSolid, box, gibbetPost, glm::vec3(0.32, 0.24, 0.15), 0.85f, 0.0f);
        r.drawPBR(r.pSolid, box, gibbetArm, glm::vec3(0.32, 0.24, 0.15), 0.85f, 0.0f);
        r.drawPBR(r.pSolid, sphere, ballM, glm::vec3(0.22, 0.21, 0.22), 0.45f, 0.95f);
        drawChain(slung, false); drawChain(flail, false); drawChain(spill, false);
        r.endScene();
    };
    auto stepAll = [&](real dt) { slung.step(dt); flail.step(dt); spill.step(dt); };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 0.35f + f * (0.45f / frames); renderFrame();
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
