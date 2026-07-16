// The Coriolis effect on a fluid, shown side by side. Two identical shallow tanks
// of water (GPU WCSPH, phys::gpu::GpuSPH) each start with a raised central mound of
// dyed water; when it collapses it drives a radial outrush. LEFT tank is in an
// inertial frame (no rotation) — the mound spreads symmetrically and its dye pinwheel
// stays radial. RIGHT tank sits in a rotating frame (a geophysical "f-plane": the
// Coriolis acceleration a = -2·omega×v is added, gravity stays vertical) — the same
// outrush is deflected sideways, spinning up a vortex that winds the pinwheel into a
// spiral. Same initial condition, one knob (omega): the swirl is purely Coriolis.
// Dye is a passive scalar advected in the particle (pos.w); drawn as a screen-space
// fluid surface from a near-top-down view.
//   ./coriolis3d [--shot out.png [f]] | --video f [n]
#include "phys/gpu/gpu_sph.cuh"
#include "common/gfx.h"
#include "common/fluidsurf.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <memory>
using namespace phys::gpu;

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 320;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    // shared fluid params (both tanks simulated in local coords centred at origin)
    const float sp = 0.009f, OMEGA = 16.0f, TANKSEP = 0.66f;
    Params base{};
    base.h = 0.020f; base.mass = 0.006f; base.stiff = 120.f; base.visc = 6.f;
    base.gmin = make_float3(-0.55f, 0.f, -0.55f); base.gmax = make_float3(0.55f, 0.8f, 0.55f);
    base.grav = make_float3(0, -9.81f, 0); base.wallDamp = 0.3f;

    // initial particles: a shallow disk of water + a raised central mound; dye is a
    // pinwheel (alternating angular sectors) so any rotation is obvious as it winds up.
    auto hash = [](float a, float b, float c) { float v = sinf(a * 127.1f + b * 311.7f + c * 74.7f) * 43758.5453f; return v - floorf(v); };
    std::vector<float3> pos; std::vector<float> dye;
    const float Rtank = 0.42f, layer = 0.09f, Rmound = 0.15f, Hmound = 0.34f;
    auto add = [&](float x, float y, float z) {
        pos.push_back(make_float3(x + (hash(x, y, z) - 0.5f) * 0.4f * sp, y + (hash(y, z, x) - 0.5f) * 0.4f * sp, z + (hash(z, x, y) - 0.5f) * 0.4f * sp));
        float ang = atan2f(z, x) + 3.14159265f;                               // [0,2pi)
        float q = fmodf(ang, 1.5707963f);                                     // position within a 90° quadrant
        dye.push_back(q < 0.72f ? 1.0f : 0.0f);                               // 4 bold dye arms (a pinwheel) that wind up under rotation
    };
    for (float x = -Rtank; x <= Rtank; x += sp) for (float z = -Rtank; z <= Rtank; z += sp) {
        float r = sqrtf(x * x + z * z); if (r > Rtank) continue;
        for (float y = 0.015f; y <= layer; y += sp) add(x, y, z);            // shallow base layer
        if (r < Rmound) { float top = layer + Hmound * (1.0f - r / Rmound);  // conical central mound
            for (float y = layer + sp; y <= top; y += sp) add(x, y, z); }
    }

    // two tanks: LEFT inertial (omega=0), RIGHT rotating (f-plane)
    Params pL = base; pL.omega = make_float3(0, 0, 0);
    Params pR = base; pR.omega = make_float3(0, OMEGA, 0);
    GpuSPH tankL, tankR; tankL.init(pos, pL, &dye); tankR.init(pos, pR, &dye);
    std::printf("coriolis: %d + %d particles (omega=%.1f)\n", tankL.n, tankR.n, OMEGA);
    const int SUB = 9; const float dt = 7.0e-4f;

    gfx::App app(W, H, "coriolis", headless);
    gfx::FluidSurf fs; fs.init(W, H);
    fs.setLight(glm::vec3(-0.35f, -0.9f, -0.25f));
    gfx::Mesh box = gfx::makeBox(), plane = gfx::makePlane(40, 1);

    gfx::OrbitCamera cam; cam.target = glm::vec3(0.0f, 0.02f, 0.0f); cam.dist = 2.5f; cam.yaw = 0.0f; cam.pitch = 1.18f;   // near-top-down
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.05f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    // tank floors (dark discs approximated by thin boxes)
    auto floorModel = [&](float cx) { return glm::scale(glm::translate(glm::mat4(1), glm::vec3(cx, -0.01f, 0)), glm::vec3(0.57f, 0.01f, 0.57f)); };

    std::vector<float4> hL, hR; std::vector<float> parts;
    auto renderFrame = [&]() {
        tankL.downloadPosDye(hL); tankR.downloadPosDye(hR);
        parts.clear(); parts.reserve((hL.size() + hR.size()) * 4);
        for (auto& q : hL) { parts.push_back(q.x - TANKSEP); parts.push_back(q.y); parts.push_back(q.z); parts.push_back(q.w); }
        for (auto& q : hR) { parts.push_back(q.x + TANKSEP); parts.push_back(q.y); parts.push_back(q.z); parts.push_back(q.w); }
        fs.setCamera(cam.view(), proj, cam.eye());
        fs.beginScene(glm::vec3(0.16, 0.17, 0.20));
        fs.drawMesh(plane, glm::mat4(1), glm::vec3(0.24, 0.24, 0.26), 0.95f, 1);
        fs.drawMesh(box, floorModel(-TANKSEP), glm::vec3(0.10, 0.10, 0.12), 0.8f, 0);
        fs.drawMesh(box, floorModel(TANKSEP), glm::vec3(0.10, 0.10, 0.12), 0.8f, 0);
        fs.surface(parts, 4, 0.019f);
        fs.blur(7, 0.05f);
        fs.composite(1, glm::vec3(0.05, 0.25, 0.60), glm::vec3(0.90, 0.15, 0.10));   // blue water ↔ red dye pinwheel
    };
    auto step = [&]() { for (int s = 0; s < SUB; s++) { tankL.step(dt); tankR.step(dt); } };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame(); char q[512]; std::snprintf(q, sizeof(q), "%s_%04d.png", video, f); fs.screenshot(q); step(); }
        std::printf("wrote %d frames\n", frames); tankL.free(); tankR.free(); return 0;
    }
    if (headless) { for (int i = 0; i < frames; i++) step(); renderFrame(); fs.screenshot(shot); std::printf("wrote %s\n", shot); tankL.free(); tankR.free(); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        step(); renderFrame(); fs.present(); app.poll();
    }
    tankL.free(); tankR.free();
    return 0;
}
