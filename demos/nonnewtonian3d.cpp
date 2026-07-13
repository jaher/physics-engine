// Non-Newtonian fluids: two open troughs get an identical heavy ball dropped in
// at the same instant. LEFT is shear-thickening "oobleck" (phys::Rheology::oobleck,
// flow index n>1), RIGHT is ordinary water. The oobleck's viscosity spikes under
// the impact's high shear rate, so it absorbs the ball almost like a solid with
// barely a splash, while the water erupts in a crown and lets the ball plunge.
// Both are drawn as smooth screen-space fluid surfaces (gfx::FluidSurf).
//   ./nonnewtonian3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "phys/sph.h"
#include "common/gfx.h"
#include "common/fluidsurf.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <memory>
using namespace phys;

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int SUB = 17;
    const real sp = 0.035, dt = 1.0e-3;

    auto makeTank = [&](real cx, Rheology rh) {
        auto f = std::make_unique<SPHFluid>();
        Vector3 cmin(cx - 0.45, 0, -0.28), cmax(cx + 0.45, 1.25, 0.28);
        f->rheo = rh; f->h = 0.072;                       // finer particles → smoother surface
        f->build(cmin, cmax, Vector3(cx - 0.4, 0.04, -0.24), Vector3(cx + 0.4, 0.62, 0.24), sp);
        return f;
    };
    std::vector<std::unique_ptr<SPHFluid>> fluids;
    fluids.push_back(makeTank(-0.62, Rheology::oobleck()));   // left: shear-thickening
    fluids.push_back(makeTank(0.62, Rheology::water()));      // right: Newtonian
    float key[2] = {0.0f, 1.0f};                             // 0 → oobleck colour, 1 → water colour

    for (auto& f : fluids) for (int s = 0; s < 500; s++) f->step(dt);          // settle the pools
    for (size_t t = 0; t < fluids.size(); t++) {                              // drop identical heavy balls
        real cx = t == 0 ? -0.62 : 0.62;
        real surf = 0; for (auto& p : fluids[t]->pos) surf = std::max(surf, (real)p.y);
        FluidBall b; b.active = true; b.radius = 0.13; b.mass = 10.0;
        b.pos = Vector3(cx, surf + 0.45, 0); b.vel = Vector3(0, -3.2, 0);
        fluids[t]->balls.push_back(b);
    }

    gfx::App app(W, H, "non-newtonian fluids", headless);
    gfx::FluidSurf fs; fs.init(W, H);
    fs.setLight(glm::vec3(-0.4f, -0.9f, -0.45f));
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 0.5, 0); cam.dist = 3.6f; cam.yaw = 0.42f; cam.pitch = 0.30f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    std::vector<glm::mat4> troughs;
    auto addTrough = [&](real cx) {
        real x0 = cx - 0.47, x1 = cx + 0.47, z0 = -0.30, z1 = 0.30, hgt = 0.68, tw = 0.03;
        troughs.push_back(glm::scale(glm::translate(glm::mat4(1), glm::vec3(cx, -0.015f, 0)), glm::vec3((x1 - x0) * 0.5f + tw, 0.02f, (z1 - z0) * 0.5f + tw)));
        troughs.push_back(glm::scale(glm::translate(glm::mat4(1), glm::vec3(cx, hgt * 0.5f, z0)), glm::vec3((x1 - x0) * 0.5f + tw, hgt * 0.5f, tw)));
        troughs.push_back(glm::scale(glm::translate(glm::mat4(1), glm::vec3(x0, hgt * 0.5f, 0)), glm::vec3(tw, hgt * 0.5f, (z1 - z0) * 0.5f)));
        troughs.push_back(glm::scale(glm::translate(glm::mat4(1), glm::vec3(x1, hgt * 0.5f, 0)), glm::vec3(tw, hgt * 0.5f, (z1 - z0) * 0.5f)));
    };
    addTrough(-0.62); addTrough(0.62);

    std::vector<float> parts;
    auto renderFrame = [&]() {
        fs.setCamera(cam.view(), proj, cam.eye());
        fs.beginScene(glm::vec3(0.46, 0.56, 0.72));
        fs.drawMesh(plane, glm::mat4(1), glm::vec3(0.34, 0.35, 0.37), 0.95f, 1);
        for (auto& m : troughs) fs.drawMesh(box, m, glm::vec3(0.11, 0.11, 0.13), 0.5f, 0);
        for (auto& f : fluids) for (auto& b : f->balls) if (b.active) {
            glm::mat4 M = glm::scale(glm::translate(glm::mat4(1), glm::vec3(b.pos.x, b.pos.y, b.pos.z)), glm::vec3((float)b.radius));
            fs.drawMesh(sphere, M, glm::vec3(0.20, 0.20, 0.23), 0.3f, 0);
        }
        parts.clear();
        for (size_t t = 0; t < fluids.size(); t++)
            for (auto& p : fluids[t]->pos) parts.insert(parts.end(), {(float)p.x, (float)p.y, (float)p.z, key[t]});
        fs.surface(parts, 4, 0.042f);
        fs.blur(5, 0.06f);
        fs.composite(1, glm::vec3(0.85, 0.82, 0.70), glm::vec3(0.10, 0.35, 0.62));   // oobleck cream, water blue
    };
    auto stepAll = [&]() { for (int s = 0; s < SUB; s++) for (auto& f : fluids) f->step(dt); };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); fs.screenshot(p); stepAll(); }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames; i++) stepAll(); renderFrame(); fs.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        stepAll(); renderFrame(); fs.present(); app.poll();
    }
    return 0;
}
