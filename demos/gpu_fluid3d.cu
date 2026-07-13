// GPU-accelerated fluid: a large-scale SPH dam-break run entirely on the GPU
// (phys::gpu::GpuSPH — CUDA grid + sorted-neighbour density/force/integrate),
// with positions read back each frame and drawn as a smooth screen-space surface.
//   ./gpu_fluid3d [--shot out.png [f]] | --video f [n]
#include "phys/gpu/gpu_sph.cuh"
#include "common/gfx.h"
#include "common/fluidsurf.h"
#include <cstring>
#include <cmath>
#include <vector>
using namespace phys::gpu;

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    // a tall column of fluid on one side of the tank → dam break.
    // Fine particles (sp≈8.5 mm, ~0.9 M of them) so the screen-space surface reads
    // as continuous liquid, not visible spheres — the GPU eats the extra count.
    Params p{}; p.h = 0.0165f; p.mass = 0.005f; p.stiff = 200.f; p.visc = 9.f;
    p.gmin = make_float3(-0.85f, 0.f, -0.4f); p.gmax = make_float3(0.85f, 1.3f, 0.4f);
    p.grav = make_float3(0, -9.81f, 0); p.wallDamp = 0.25f;
    const float sp = 0.0085f;
    auto hash = [](float a, float b, float c) { float v = sinf(a * 127.1f + b * 311.7f + c * 74.7f) * 43758.5453f; return v - floorf(v); };
    std::vector<float3> pos;
    for (float x = -0.82f; x <= -0.1f; x += sp) for (float y = 0.03f; y <= 1.05f; y += sp) for (float z = -0.37f; z <= 0.37f; z += sp)
        // jitter off the perfect lattice so the flat surface doesn't shimmer with grid moiré
        pos.push_back(make_float3(x + (hash(x, y, z) - 0.5f) * 0.4f * sp,
                                  y + (hash(y, z, x) - 0.5f) * 0.4f * sp,
                                  z + (hash(z, x, y) - 0.5f) * 0.4f * sp));
    GpuSPH sph; sph.init(pos, p);
    std::printf("GPU dam-break: %d particles\n", sph.n);

    gfx::App app(W, H, "gpu fluid", headless);
    gfx::FluidSurf fs; fs.init(W, H);
    fs.setLight(glm::vec3(-0.4f, -0.9f, -0.45f));
    gfx::Mesh box = gfx::makeBox();
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 0.35, 0); cam.dist = 3.0f; cam.yaw = 0.5f; cam.pitch = 0.24f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    float hx = 0.85f, hz = 0.4f, hy = 1.3f, tw = 0.02f;
    glm::mat4 back = glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, hy * 0.5f, -hz - tw)), glm::vec3(hx + tw, hy * 0.5f, tw));
    glm::mat4 lft = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-hx - tw, hy * 0.5f, 0)), glm::vec3(tw, hy * 0.5f, hz));
    glm::mat4 rgt = glm::scale(glm::translate(glm::mat4(1), glm::vec3(hx + tw, hy * 0.5f, 0)), glm::vec3(tw, hy * 0.5f, hz));

    std::vector<float3> host; std::vector<float> parts;
    auto renderFrame = [&]() {
        sph.download(host);
        parts.clear(); parts.reserve(host.size() * 4);
        for (auto& q : host) { parts.push_back(q.x); parts.push_back(q.y); parts.push_back(q.z); parts.push_back(0); }
        fs.setCamera(cam.view(), proj, cam.eye());
        fs.beginScene(glm::vec3(0.46, 0.56, 0.72));
        fs.drawMesh(box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, -0.02f, 0)), glm::vec3(hx + tw, 0.02f, hz + tw)), glm::vec3(0.25, 0.26, 0.30), 0.9f, 1);
        fs.drawMesh(box, back, glm::vec3(0.16, 0.17, 0.2), 0.6f, 0);
        fs.drawMesh(box, lft, glm::vec3(0.16, 0.17, 0.2), 0.6f, 0);
        fs.drawMesh(box, rgt, glm::vec3(0.16, 0.17, 0.2), 0.6f, 0);
        fs.surface(parts, 4, 0.017f);      // small imposters (≈2·spacing) → tiny bumps
        fs.blur(8, 0.06f);                 // bilateral smoothing fuses them into one surface
        fs.composite(1, glm::vec3(0.12, 0.4, 0.7));
    };
    auto step = [&]() { for (int s = 0; s < 16; s++) sph.step(4.0e-4f); };   // smaller h ⇒ smaller stable dt

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame(); char q[512]; std::snprintf(q, sizeof(q), "%s_%04d.png", video, f); fs.screenshot(q); step(); }
        std::printf("wrote %d frames\n", frames); sph.free(); return 0;
    }
    if (headless) { for (int i = 0; i < frames; i++) step(); renderFrame(); fs.screenshot(shot); std::printf("wrote %s\n", shot); sph.free(); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        step(); renderFrame(); fs.present(); app.poll();
    }
    sph.free();
    return 0;
}
