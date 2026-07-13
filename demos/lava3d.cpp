// Hyper-real lava flowing around a cylinder. A hot, shear-thinning Bingham fluid
// (phys::Rheology::lava — high viscosity + yield stress) is poured onto a tilted
// channel (the slope is emulated with a down-slope gravity component) and streams
// past a vertical stone cylinder, splitting and rejoining in its wake. Each SPH
// particle carries a temperature that cools with age; the fluid is drawn as a
// smooth screen-space surface (gfx::FluidSurf) with a molten blackbody emissive —
// bright yellow-white at the vent, darkening to a red-black crust downstream.
//   ./lava3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "phys/sph.h"
#include "common/gfx.h"
#include "common/fluidsurf.h"
#include <cstring>
#include <cmath>
#include <vector>
using namespace phys;

static gfx::Mesh makeCylinder(int seg = 64) {
    std::vector<float> v; std::vector<unsigned> idx;
    auto push = [&](float x, float y, float z, float nx, float ny, float nz) { v.insert(v.end(), {x, y, z, nx, ny, nz}); };
    for (int i = 0; i <= seg; i++) { float a = 2 * (float)M_PI * i / seg, c = std::cos(a), s = std::sin(a); push(c, 0, s, c, 0, s); push(c, 1, s, c, 0, s); }
    for (int i = 0; i < seg; i++) { unsigned b = i * 2; idx.insert(idx.end(), {b, b + 1, b + 2, b + 2, b + 1, b + 3}); }
    unsigned center = (unsigned)(v.size() / 6); push(0, 1, 0, 0, 1, 0);
    unsigned ring0 = (unsigned)(v.size() / 6);
    for (int i = 0; i <= seg; i++) { float a = 2 * (float)M_PI * i / seg; push(std::cos(a), 1, std::sin(a), 0, 1, 0); }
    for (int i = 0; i < seg; i++) idx.insert(idx.end(), {center, ring0 + (unsigned)i, ring0 + (unsigned)i + 1});
    gfx::Mesh m; m.upload(v, idx); return m;
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int SUB = 17;
    const real dt = 1.0e-3, sp = 0.030;
    const real slope = 30.0 * 3.14159265 / 180.0;

    SPHFluid lava;
    Vector3 cmin(-1.25, 0, -0.34), cmax(1.25, 0.7, 0.34);
    lava.rheo = Rheology::lava();
    lava.h = 0.062;                                       // finer particles → smoother surface
    lava.gravity = Vector3(9.81 * std::sin(slope), -9.81 * std::cos(slope), 0);
    lava.restitution = 0.05; lava.wallFriction = 0.5; lava.xsphEps = 0.22;
    lava.build(cmin, cmax, Vector3(-1.18, 0.04, -0.28), Vector3(-0.72, 0.4, 0.28), sp);
    SPHCylinder cyl; cyl.cx = 0.15; cyl.cz = 0.0; cyl.radius = 0.17; cyl.y0 = 0; cyl.y1 = 1e9; cyl.friction = 0.6;
    lava.cylinders.push_back(cyl);

    std::vector<float> temp(lava.nParticles(), 1.0f);
    unsigned seed = 99; auto rnd = [&]() { seed = seed * 1103515245u + 12345u; return (float)(((seed >> 16) & 0x7fff) / 32767.0); };
    auto emit = [&]() {
        for (real z = -0.26; z <= 0.26; z += sp)
            for (real y = 0.26; y <= 0.40; y += sp) {
                lava.addParticle(Vector3(-1.18 + 0.02 * rnd(), y, z + 0.01 * rnd()), Vector3(1.6, -0.2, 0));
                temp.push_back(1.0f);
            }
    };
    auto drain = [&]() {
        for (int i = 0; i < lava.nParticles();) {
            if (lava.pos[i].x > cmax.x - 0.05) { lava.removeParticle(i); temp[i] = temp.back(); temp.pop_back(); }
            else i++;
        }
    };

    gfx::App app(W, H, "lava", headless);
    gfx::FluidSurf fs; fs.init(W, H);
    fs.setLight(glm::vec3(-0.35f, -0.85f, -0.4f));
    gfx::Mesh plane = gfx::makePlane(40, 1), box = gfx::makeBox(), cylMesh = makeCylinder();

    gfx::OrbitCamera cam; cam.target = glm::vec3(0.0, 0.05, 0); cam.dist = 3.05f; cam.yaw = 0.62f; cam.pitch = 0.46f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    glm::mat4 cylM = glm::scale(glm::translate(glm::mat4(1), glm::vec3((float)cyl.cx, 0.0f, (float)cyl.cz)), glm::vec3((float)cyl.radius, 0.62f, (float)cyl.radius));
    glm::mat4 floorM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, -0.02f, 0)), glm::vec3(1.3f, 0.02f, 0.36f));

    std::vector<float> parts;
    auto renderFrame = [&]() {
        fs.setCamera(cam.view(), proj, cam.eye());
        fs.beginScene(glm::vec3(0.04, 0.035, 0.045));
        fs.drawMesh(plane, glm::mat4(1), glm::vec3(0.05, 0.045, 0.05), 0.95f, 1);
        fs.drawMesh(box, floorM, glm::vec3(0.08, 0.07, 0.07), 0.9f, 0);
        fs.drawMesh(cylMesh, cylM, glm::vec3(0.14, 0.13, 0.13), 0.7f, 0);
        parts.clear();
        for (int i = 0; i < lava.nParticles(); i++) { const Vector3& p = lava.pos[i]; parts.insert(parts.end(), {(float)p.x, (float)p.y, (float)p.z, temp[i]}); }
        fs.surface(parts, 4, 0.036f);
        fs.blur(5, 0.05f);
        fs.composite(0, glm::vec3(0));
    };
    auto stepAll = [&]() {
        emit();
        for (int s = 0; s < SUB; s++) lava.step(dt);
        drain();
        float cool = std::exp(-0.85f * (SUB * (float)dt));
        for (auto& t : temp) t *= cool;
    };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); fs.screenshot(p); stepAll(); }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames; i++) stepAll(); renderFrame(); fs.screenshot(shot); std::printf("wrote %s (N=%d)\n", shot, lava.nParticles()); return 0; }
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
