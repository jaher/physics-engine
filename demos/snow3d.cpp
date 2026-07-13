// Snow with footprints: a walking robot crosses a fresh snowfield; each footfall
// presses a print into the deformable heightfield (HeightField::deform), leaving
// a trail behind it. Simple kinematic biped gait; snow shaded bright and soft.
//   ./snow3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
using namespace phys;

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int GN = 200;                                    // fine grid so prints are crisp

    HeightField snow; snow.init(GN, GN, 0.08, Vector3(-GN * 0.04, 0, -GN * 0.04));
    unsigned rs = 42; auto rr = [&]() { rs = rs * 1103515245u + 12345u; return ((rs >> 16) & 0x7fff) / 32767.0; };
    for (int z = 0; z < GN; z++) for (int x = 0; x < GN; x++)
        snow.at(x, z) = 0.03 * std::sin(x * 0.2) * std::cos(z * 0.18) + 0.008 * rr();   // windswept sparkle

    // ---- kinematic biped gait along a curved path
    auto pathPos = [&](double s) { return Vector3(-6.5 + s, 0, -2.0 + 1.5 * std::sin(s * 0.45)); };
    const real stepLen = 0.55, stepTime = 0.42, footR = 0.16, hipW = 0.20;
    double walkT = 0; int stepIndex = 0; bool leftSwing = true;
    Vector3 plantedL, plantedR;                            // world positions of planted feet
    {
        Vector3 p0 = pathPos(0); Vector3 fwd = (pathPos(0.01) - p0).unit(); Vector3 side = Vector3(0, 1, 0) % fwd;
        plantedL = p0 + side * hipW; plantedR = p0 - side * hipW;
        plantedL.y = snow.sample(plantedL.x, plantedL.z); plantedR.y = snow.sample(plantedR.x, plantedR.z);
    }
    Vector3 swingFrom = plantedL, swingTo = plantedL;
    bool stampedThisStep = true;

    gfx::App app(W, H, "snow footprints", headless);
    gfx::Renderer r; r.init(W, H);
    r.lightColor = glm::vec3(1.0f, 0.97f, 0.9f) * 3.0f;
    r.lightDir = glm::normalize(glm::vec3(-0.95f, -0.42f, -0.25f));   // low winter sun rakes the relief
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere();
    gfx::Mesh snowMesh;
    std::vector<float> tv; std::vector<unsigned> ti;
    auto buildSnow = [&](bool first) {
        tv.clear();
        for (int z = 0; z < GN; z++) for (int x = 0; x < GN; x++) {
            real wx = snow.origin.x + x * snow.spacing, wz = snow.origin.z + z * snow.spacing;
            Vector3 n = snow.normal(wx, wz);
            tv.insert(tv.end(), {(float)wx, (float)snow.sample(wx, wz), (float)wz, (float)n.x, (float)n.y, (float)n.z});
        }
        if (first) { ti.clear();
            for (int z = 0; z < GN - 1; z++) for (int x = 0; x < GN - 1; x++) {
                unsigned a = z * GN + x, b = a + 1, c = a + GN, d = c + 1;
                ti.insert(ti.end(), {a, c, b, b, c, d}); }
            snowMesh.upload(tv, ti); }
        else snowMesh.updateVerts(tv);
    };
    buildSnow(true);

    gfx::OrbitCamera cam; cam.target = glm::vec3(-1.5, 0.7, -1.5); cam.dist = 6.2f; cam.yaw = 2.45f; cam.pitch = 0.5f;
    glm::mat4 proj = glm::perspective(glm::radians(43.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    Vector3 bodyPos, footL, footR_;
    auto stepGait = [&](double dt) {
        walkT += dt;
        double phase = std::fmod(walkT, (double)stepTime) / stepTime;
        double s = walkT / stepTime * stepLen * 0.5;       // body travel
        Vector3 p = pathPos(s); Vector3 fwd = (pathPos(s + 0.01) - p).unit(); Vector3 side = Vector3(0, 1, 0) % fwd;
        int idx = (int)(walkT / stepTime);
        if (idx != stepIndex) {                            // start a new step
            stepIndex = idx; leftSwing = !leftSwing; stampedThisStep = false;
            Vector3& sw = leftSwing ? plantedL : plantedR;
            swingFrom = sw;
            double sAhead = s + stepLen * 1.6;
            Vector3 ahead = pathPos(sAhead); Vector3 f2 = (pathPos(sAhead + 0.01) - ahead).unit(); Vector3 sd2 = Vector3(0, 1, 0) % f2;
            swingTo = ahead + sd2 * (leftSwing ? hipW : -hipW);
            swingTo.y = snow.sample(swingTo.x, swingTo.z);
        }
        // swing foot arcs from→to; stance foot stays planted
        Vector3 swing = swingFrom + (swingTo - swingFrom) * phase;
        swing.y = swingFrom.y + (swingTo.y - swingFrom.y) * phase + 0.28 * std::sin(phase * M_PI);
        if (leftSwing) { footL = swing; footR_ = plantedR; } else { footR_ = swing; footL = plantedL; }
        if (phase > 0.92 && !stampedThisStep) {            // footfall → stamp the print
            Vector3& sw = leftSwing ? plantedL : plantedR;
            sw = swingTo;
            snow.deform(swingTo.x, swingTo.z, footR * 2.0, footR * 3.0, 0.13);
            stampedThisStep = true; buildSnow(false);
        }
        Vector3 mid = (footL + footR_) * 0.5;
        bodyPos = Vector3(mid.x, snow.sample(p.x, p.z) + 1.05 + 0.03 * std::sin(walkT * 7), mid.z);
    };

    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 0.5, 0), 10.0f);
        Vector3 fwd = (pathPos(walkT / stepTime * stepLen * 0.5 + 0.01) - pathPos(walkT / stepTime * stepLen * 0.5)).unit();
        float yawA = (float)std::atan2(fwd.x, fwd.z);
        auto place = [&](Vector3 p, glm::vec3 sc) { return glm::scale(glm::rotate(glm::translate(glm::mat4(1), glm::vec3(p.x, p.y, p.z)), yawA, glm::vec3(0, 1, 0)), sc); };
        glm::mat4 torso = place(bodyPos, glm::vec3(0.28f, 0.42f, 0.18f));
        glm::mat4 head = place(bodyPos + Vector3(0, 0.62, 0), glm::vec3(0.16f));
        glm::mat4 fl = place(footL + Vector3(0, footR * 0.5, 0), glm::vec3(0.11f, 0.07f, 0.19f));
        glm::mat4 fr2 = place(footR_ + Vector3(0, footR * 0.5, 0), glm::vec3(0.11f, 0.07f, 0.19f));
        auto leg = [&](Vector3 foot, real sideSign) {
            Vector3 hip = bodyPos + Vector3(0, -0.42, 0);
            Vector3 mid2 = (hip + foot) * 0.5; Vector3 d = hip - foot;
            glm::mat4 m = glm::translate(glm::mat4(1), glm::vec3(mid2.x, mid2.y, mid2.z));
            return glm::scale(m, glm::vec3(0.06f, (float)(d.magnitude() * 0.5), 0.06f)); (void)sideSign; };
        r.beginShadow();
        r.shadowDraw(snowMesh, glm::mat4(1));
        r.shadowDraw(box, torso); r.shadowDraw(sphere, head); r.shadowDraw(box, fl); r.shadowDraw(box, fr2);
        r.shadowDraw(box, leg(footL, 1)); r.shadowDraw(box, leg(footR_, -1));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pSolid, snowMesh, glm::mat4(1), glm::vec3(0.93, 0.95, 0.99), 0.55f, 0.0f);   // bright soft snow
        r.drawPBR(r.pSolid, box, torso, glm::vec3(0.75, 0.30, 0.12), 0.45f, 0.4f);
        r.drawPBR(r.pSolid, sphere, head, glm::vec3(0.85, 0.86, 0.9), 0.3f, 0.7f);
        r.drawPBR(r.pSolid, box, fl, glm::vec3(0.2, 0.2, 0.24), 0.6f, 0.3f);
        r.drawPBR(r.pSolid, box, fr2, glm::vec3(0.2, 0.2, 0.24), 0.6f, 0.3f);
        r.drawPBR(r.pSolid, box, leg(footL, 1), glm::vec3(0.5, 0.52, 0.56), 0.5f, 0.5f);
        r.drawPBR(r.pSolid, box, leg(footR_, -1), glm::vec3(0.5, 0.52, 0.56), 0.5f, 0.5f);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 2.45f - f * (0.4f / frames);
            renderFrame(); char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 2; s++) stepGait(1.0 / 120); }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames * 2; i++) stepGait(1.0 / 120); renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        stepGait(1.0 / 60); renderFrame(); r.present(); app.poll();
    }
    return 0;
}
