// Robotics layer (phys::robotics): a serial arm tracks a moving target with the
// Jacobian-transpose IK solver, while a Lidar sensor sweeps a fan of rays over the
// obstacles around its base. Reduced-coordinate forward kinematics places every
// link; the same header also has IMU/contact-force sensors and RNE inverse dynamics.
//   ./robotics3d [--shot out.png [f]] | --video f [n]
#include "phys/phys.h"
#include "phys/robotics.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace phys;

// place a unit box ([-1,1]³, long axis = local Y) as a rod from a to b
static glm::mat4 rod(glm::vec3 a, glm::vec3 b, float rad) {
    glm::vec3 mid = 0.5f * (a + b), d = b - a; float len = glm::length(d);
    glm::vec3 dir = len > 1e-5f ? d / len : glm::vec3(0, 1, 0);
    glm::vec3 up = std::fabs(dir.y) < 0.98f ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
    glm::vec3 rt = glm::normalize(glm::cross(up, dir)), fw = glm::cross(dir, rt);
    glm::mat4 R(glm::vec4(rt, 0), glm::vec4(dir, 0), glm::vec4(fw, 0), glm::vec4(0, 0, 0, 1));
    return glm::translate(glm::mat4(1), mid) * R * glm::scale(glm::mat4(1), glm::vec3(rad, len * 0.5f, rad));
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    SerialChain chain; chain.baseOrigin = Vector3(-1.0, 0.7, 0);
    for (int i = 0; i < 4; i++) { ChainLink L; L.axis = Vector3(0, 0, 1); L.toChild = Vector3(0.52, 0, 0); L.com = Vector3(0.26, 0, 0); chain.links.push_back(L); }
    std::vector<real> q = {0.6, -0.5, -0.5, -0.4};
    JacobianIK ik; ik.chain = &chain;

    // lidar obstacles: vertical box pillars (AABB footprint in x-z)
    struct Pillar { real x, z, hx, hz; }; std::vector<Pillar> pillars = {{0.7, -0.9, 0.12, 0.12}, {1.4, 0.5, 0.14, 0.14}, {-0.3, 1.1, 0.1, 0.1}, {0.2, -1.5, 0.13, 0.13}};
    Lidar lidar; lidar.numRays = 96; lidar.fov = 6.2831853; lidar.maxRange = 3.2; lidar.origin = Vector3(-1.0, 0.14, 0); lidar.forward = Vector3(1, 0, 0); lidar.up = Vector3(0, 1, 0);
    auto raycast = [&](Vector3 o, Vector3 dir) -> real {                  // 2D ray vs AABB pillars
        real best = lidar.maxRange;
        for (auto& p : pillars) {
            real tmin = 0, tmax = best; bool hit = true;
            real bo[2] = {o.x, o.z}, bd[2] = {dir.x, dir.z}, lo[2] = {p.x - p.hx, p.z - p.hz}, hi[2] = {p.x + p.hx, p.z + p.hz};
            for (int a = 0; a < 2; a++) { if (std::fabs(bd[a]) < 1e-8) { if (bo[a] < lo[a] || bo[a] > hi[a]) { hit = false; break; } }
                else { real t1 = (lo[a] - bo[a]) / bd[a], t2 = (hi[a] - bo[a]) / bd[a]; if (t1 > t2) std::swap(t1, t2); tmin = std::max(tmin, t1); tmax = std::min(tmax, t2); if (tmin > tmax) { hit = false; break; } } }
            if (hit && tmin >= 0) best = std::min(best, tmin);
        }
        return best;
    };

    gfx::App app(W, H, "robot arm + lidar", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0.1, 1.1, 0); cam.dist = 5.6f; cam.yaw = 0.5f; cam.pitch = 0.26f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    double t = 0; std::vector<real> ranges;
    auto step = [&]() {
        t += 1.0 / 60;
        Vector3 target(0.35 + 0.75 * std::cos(t * 0.9), 1.5 + 0.55 * std::sin(t * 1.3), 0);
        for (int k = 0; k < 8; k++) ik.step(q, target);                  // chase the target
        ranges = lidar.scan(raycast);
    };
    Vector3 target(1.1, 1.5, 0);

    auto renderFrame = [&]() {
        target = Vector3(0.35 + 0.75 * std::cos(t * 0.9), 1.5 + 0.55 * std::sin(t * 1.3), 0);
        FKResult fk = chain.fk(q);
        r.setLightForScene(glm::vec3(0, 1.1, 0), 5.5f);
        r.beginShadow();
        for (int i = 0; i < (int)chain.links.size(); i++) { glm::vec3 a(fk.jointPos[i].x, fk.jointPos[i].y, fk.jointPos[i].z), b(fk.jointPos[i + 1].x, fk.jointPos[i + 1].y, fk.jointPos[i + 1].z);
            r.shadowDraw(box, rod(a, b, 0.05f)); }
        for (auto& p : pillars) r.shadowDraw(box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(p.x, 0.5f, p.z)), glm::vec3(p.hx, 0.5f, p.hz)));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.20, 0.21, 0.24), 0.95f, 0.0f);
        // pedestal + arm links
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(-1.0f, 0.35f, 0)), glm::vec3(0.12f, 0.35f, 0.12f)), glm::vec3(0.15, 0.15, 0.17), 0.6f, 0.3f);
        for (int i = 0; i < (int)chain.links.size(); i++) {
            glm::vec3 a(fk.jointPos[i].x, fk.jointPos[i].y, fk.jointPos[i].z), b(fk.jointPos[i + 1].x, fk.jointPos[i + 1].y, fk.jointPos[i + 1].z);
            float u = (float)i / 3; glm::vec3 col = glm::mix(glm::vec3(0.85, 0.5, 0.15), glm::vec3(0.2, 0.5, 0.75), u);
            r.drawPBR(r.pSolid, box, rod(a, b, 0.05f), col, 0.35f, 0.4f);
            r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), a), glm::vec3(0.075f)), glm::vec3(0.3, 0.31, 0.34), 0.3f, 0.8f);
        }
        // end-effector + target
        glm::vec3 ee(fk.endEffector.x, fk.endEffector.y, fk.endEffector.z), tg(target.x, target.y, target.z);
        r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), ee), glm::vec3(0.08f)), glm::vec3(0.9, 0.85, 0.3), 0.25f, 0.6f);
        r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), tg), glm::vec3(0.09f)), glm::vec3(0.9, 0.25, 0.3), 0.4f, 0.0f);
        // pillars
        for (auto& p : pillars) r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(p.x, 0.5f, p.z)), glm::vec3(p.hx, 0.5f, p.hz)), glm::vec3(0.28, 0.26, 0.24), 0.7f, 0.1f);
        // lidar rays
        if (!ranges.empty()) {
            Vector3 f = lidar.forward.unit(), left = (lidar.up.unit() % f).unit();
            for (int i = 0; i < (int)ranges.size(); i++) {
                real a = -lidar.fov * 0.5 + lidar.fov * i / (ranges.size() - 1);
                Vector3 dir = (f * std::cos(a) + left * std::sin(a)).unit();
                Vector3 hit = lidar.origin + dir * ranges[i];
                glm::vec3 o(lidar.origin.x, lidar.origin.y, lidar.origin.z), hp(hit.x, hit.y, hit.z);
                bool near = ranges[i] < lidar.maxRange - 1e-3;
                r.drawPBR(r.pSolid, box, rod(o, hp, 0.005f), near ? glm::vec3(0.3, 0.9, 0.55) : glm::vec3(0.2, 0.4, 0.5), 0.4f, 0.0f);
                if (near) r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), hp), glm::vec3(0.02f)), glm::vec3(0.4, 1.0, 0.6), 0.3f, 0.0f);
            }
        }
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame(); char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p); step(); }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames; i++) step(); renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        step(); renderFrame(); r.present(); app.poll();
    }
    return 0;
}
