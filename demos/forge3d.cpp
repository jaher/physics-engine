// Forging demo: an iron bar clamped in a vice is (1) pressed cold — it springs
// back; (2) heated by a torch until it glows; (3) pressed again — the hot metal
// yields and takes a permanent bend; (4) the torch is removed and the bar cools,
// keeping its new shape (phys::PlasticRod, elastoplastic + thermal softening).
// Heat is rendered as incandescent glow (black→red→orange→yellow, HDR).
//   ./forge3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
using namespace phys;

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 340;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int N = 26; const real SEG = 0.105, ROD_R = 0.045;

    PlasticRod rod; rod.build(Vector3(-1.35, 1.5, 0), Vector3(1, 0, 0), N, SEG);
    rod.pin(0); rod.pin(1); rod.pin(2);

    // timeline (seconds): cold press → heat → hot press → cool/release
    auto phase = [](double t) {
        if (t < 2.2) return 0;          // cold press
        if (t < 3.0) return 1;          // release (springback)
        if (t < 7.0) return 2;          // torch on
        if (t < 10.5) return 3;         // hot press (bends!)
        return 4;                       // release + cool: permanent set
    };

    gfx::App app(W, H, "forging iron", headless);
    gfx::Renderer r; r.init(W, H);
    r.lightDir = glm::normalize(glm::vec3(-0.5f, -0.9f, -0.45f));
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(30, 1);

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.25, 0); cam.dist = 4.6f; cam.yaw = 0.35f; cam.pitch = 0.22f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    double t = 0;
    auto stepPhysics = [&](double dt) {
        int ph = phase(t);
        Vector3 tool;
        tool = ph == 0 ? Vector3(0, -160, 0) : (ph == 3 ? Vector3(0, -320, 0) : Vector3());   // cold press gentler (stays elastic)
        if (ph == 2 || ph == 3) rod.heat(N * 0.55, 6.0, 1.6, dt);    // torch on through the bend
        rod.step(dt, tool, N);
        t += dt;
    };
    // incandescent colour ramp (HDR values >1 glow through the ACES tonemap)
    auto glowColour = [&](real tp) {
        glm::vec3 iron(0.32f, 0.33f, 0.36f);
        glm::vec3 red(1.6f, 0.15f, 0.05f), orange(2.6f, 1.0f, 0.12f), yellow(3.4f, 2.6f, 0.9f);
        float x = (float)tp;
        glm::vec3 g = x < 0.45f ? glm::mix(iron, red, x / 0.45f)
                    : (x < 0.8f ? glm::mix(red, orange, (x - 0.45f) / 0.35f)
                                : glm::mix(orange, yellow, (x - 0.8f) / 0.2f));
        return g;
    };
    auto segMat = [&](int i) {                                        // oriented box for segment i
        Vector3 a = rod.pos[i], b = rod.pos[i + 1];
        Vector3 mid = (a + b) * 0.5, d = b - a;
        float len = (float)d.magnitude(); d.normalise();
        glm::vec3 dir((float)d.x, (float)d.y, (float)d.z);
        glm::vec3 z(0, 0, 1); glm::vec3 axis = glm::cross(z, dir);
        float ang = std::acos(glm::clamp(glm::dot(z, dir), -1.0f, 1.0f));
        glm::mat4 m = glm::translate(glm::mat4(1), glm::vec3(mid.x, mid.y, mid.z));
        if (glm::length(axis) > 1e-6f) m = glm::rotate(m, ang, glm::normalize(axis));
        return glm::scale(m, glm::vec3(ROD_R, ROD_R, len * 0.55f));
    };
    auto renderFrame = [&]() {
        int ph = phase(t);
        r.setLightForScene(glm::vec3(0, 1.2, 0), 4.5f);
        // vice + pedestal + floor
        glm::mat4 viceA = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-1.5f, 1.5f, 0)), glm::vec3(0.22f, 0.3f, 0.2f));
        glm::mat4 viceB = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-1.5f, 0.75f, 0)), glm::vec3(0.12f, 0.75f, 0.12f));
        r.beginShadow();
        r.shadowDraw(box, viceA); r.shadowDraw(box, viceB);
        for (int i = 0; i < N; i++) r.shadowDraw(box, segMat(i));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.35, 0.33, 0.31), 0.9f, 0.0f);
        r.drawPBR(r.pSolid, box, viceA, glm::vec3(0.2, 0.21, 0.24), 0.6f, 0.6f);
        r.drawPBR(r.pSolid, box, viceB, glm::vec3(0.2, 0.21, 0.24), 0.6f, 0.6f);
        for (int i = 0; i < N; i++) {                                 // the bar, glowing where hot
            real tp = (rod.temp[i] + rod.temp[i + 1]) * 0.5;
            r.drawPBR(r.pSolid, box, segMat(i), glowColour(tp), 0.45f - 0.25f * (float)tp, 0.85f * (1 - (float)tp));
        }
        if (ph == 2 || ph == 3) {                                     // torch flame under the hot spot
            Vector3 hp = rod.pos[(int)(N * 0.55)];
            for (int k = 0; k < 3; k++) {
                float s = 0.10f - k * 0.025f;
                glm::vec3 fc = k == 0 ? glm::vec3(0.3f, 0.6f, 2.8f) : (k == 1 ? glm::vec3(1.2f, 1.6f, 3.2f) : glm::vec3(3.0f, 3.0f, 3.4f));
                r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1),
                    glm::vec3(hp.x, hp.y - 0.18f - 0.05f * k, hp.z)), glm::vec3(s, s * 1.7f, s)), fc, 0.9f, 0.0f);
            }
            glm::mat4 nozzle = glm::scale(glm::translate(glm::mat4(1), glm::vec3(hp.x, hp.y - 0.55f, hp.z)), glm::vec3(0.045f, 0.16f, 0.045f));
            r.drawPBR(r.pSolid, box, nozzle, glm::vec3(0.15, 0.16, 0.18), 0.4f, 0.8f);
        }
        if (ph == 0 || ph == 3) {                                     // press tool above the tip
            Vector3 tip = rod.pos[N];
            glm::mat4 tool = glm::scale(glm::translate(glm::mat4(1), glm::vec3(tip.x, tip.y + 0.16f, tip.z)), glm::vec3(0.09f, 0.12f, 0.09f));
            r.drawPBR(r.pSolid, box, tool, glm::vec3(0.45, 0.3, 0.15), 0.7f, 0.1f);
        }
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 8; s++) stepPhysics(1.0 / 200); }     // 13.6 s of sim over the video
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames * 8; i++) stepPhysics(1.0 / 200); renderFrame(); r.screenshot(shot); std::printf("wrote %s (t=%.1fs set=%.2f rad)\n", shot, t, (double)rod.permanentSetAngle()); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        for (int s = 0; s < 3; s++) stepPhysics(1.0 / 200);
        renderFrame(); r.present(); app.poll();
    }
    return 0;
}
