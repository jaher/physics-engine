// A waterfall with real fluid dynamics: a GPU WCSPH fluid (phys::gpu::GpuSPH —
// CUDA uniform-grid density → pressure → viscosity → surface-tension → integrate,
// ~2.4M particles) pours off a clifftop reservoir and plunges into a pool that
// overflows its front lip. Water is a *recirculating* source/sink: particles that
// drain past the front spill respawn across a wide clifftop reservoir, so a fixed
// particle budget sustains a continuous fall. The cliff, pool floor and retaining
// lip are solid box (SDF) obstacles the fluid collides with; surface-tension
// cohesion holds the falling curtain together. Rendered as a screen-space liquid
// surface (gfx::FluidSurf) for the dense/calm water (blue pool + reservoir), plus a
// whitewater foam-mist billboard pass for the fast, aerated water (split per-particle
// by density + speed) — solid blue at top, misty white curtain, blue plunge pool.
//   ./waterfall3d [--shot out.png [f]] | --video f [n]
#include "phys/gpu/gpu_sph.cuh"
#include "common/gfx.h"
#include "common/fluidsurf.h"
#include <cstring>
#include <cmath>
#include <vector>
using namespace phys::gpu;

// axis-aligned solid box: [mn,mx] world AABB, drawn as rock, collided by the fluid.
struct Box { float3 mn, mx; };
static glm::mat4 boxModel(const Box& b) {
    glm::vec3 c(0.5f * (b.mn.x + b.mx.x), 0.5f * (b.mn.y + b.mx.y), 0.5f * (b.mn.z + b.mx.z));
    glm::vec3 h(0.5f * (b.mx.x - b.mn.x), 0.5f * (b.mx.y - b.mn.y), 0.5f * (b.mx.z - b.mn.z));
    return glm::scale(glm::translate(glm::mat4(1), c), h);
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 420;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    // ---- fluid parameters ----
    Params p{};
    p.h = 0.016f; p.mass = 0.004f; p.stiff = 180.f; p.visc = 6.f;      // low viscosity → lively, splashy water
    p.surfTens = 25.f;                                                  // just a touch of cohesion — fine splash droplets, not clumps or syrup
    p.gmin = make_float3(-0.95f, -0.5f, -1.0f); p.gmax = make_float3(0.95f, 2.9f, 1.0f);
    p.grav = make_float3(0, -9.81f, 0); p.wallDamp = 0.2f;

    // ---- solid rock geometry (clifftop channel → free drop → plunge pool) ----
    std::vector<Box> boxes = {
        {{-0.95f, -0.5f, -1.00f}, {0.95f, 2.40f, -0.35f}},   // cliff: tall back block, front face z=-0.35, top y=2.40
        {{-0.95f,  0.00f, -0.35f}, {0.95f, 0.16f, 0.78f}},   // pool floor: top y=0.16, front edge z=0.78
        {{-0.95f,  0.16f,  0.72f}, {0.95f, 0.72f, 0.80f}},   // front lip: deep pool retaining wall, top y=0.72
    };
    p.nBox = (int)boxes.size(); p.boxDamp = 0.2f; p.boxFric = 0.6f;
    for (int i = 0; i < p.nBox; i++) { p.boxMin[i] = boxes[i].mn; p.boxMax[i] = boxes[i].mx; }

    // ---- recirculating emitter: a WIDE shallow reservoir on the clifftop that
    //      overflows the front lip. Respawning drained particles across this broad,
    //      low-density channel (not a narrow slot) buffers the bursty recycle into a
    //      smooth pour instead of a geyser ----
    p.emitOn = 1;
    p.emitMin = make_float3(-0.82f, 2.40f, -0.95f); p.emitMax = make_float3(0.82f, 2.55f, -0.45f);
    p.emitVel = make_float3(0.f, 0.f, 0.35f);           // gentle push over the lip so the curtain arcs off the rock face
    p.recycleY = -0.20f;                                // drained past the spill → respawn upstream

    // ---- initial particles: clifftop reservoir + plunge pool + a hanging sheet ----
    const float sp = 0.008f;    // fine spacing → ~2.4M particles → smooth surface, no visible blobs
    auto hash = [](float a, float b, float c) { float v = sinf(a * 127.1f + b * 311.7f + c * 74.7f) * 43758.5453f; return v - floorf(v); };
    std::vector<float3> pos;
    auto fill = [&](float3 lo, float3 hi) {
        for (float x = lo.x; x <= hi.x; x += sp) for (float y = lo.y; y <= hi.y; y += sp) for (float z = lo.z; z <= hi.z; z += sp)
            pos.push_back(make_float3(x + (hash(x, y, z) - 0.5f) * 0.4f * sp, y + (hash(y, z, x) - 0.5f) * 0.4f * sp, z + (hash(z, x, y) - 0.5f) * 0.4f * sp));
    };
    fill(make_float3(-0.82f, 2.42f, -0.93f), make_float3(0.82f, 2.66f, -0.42f));  // clifftop feed reservoir
    fill(make_float3(-0.88f, 0.17f, -0.30f), make_float3(0.88f, 0.62f, 0.68f));   // deep plunge pool (below the lip)
    fill(make_float3(-0.48f, 0.30f, -0.33f), make_float3(0.48f, 2.30f, -0.22f));  // sheet seeding the initial fall

    GpuSPH sph; sph.init(pos, p);
    std::printf("waterfall: %d particles\n", sph.n);
    const int SUB = 42; const float dt = 1.8e-4f;
    for (int s = 0; s < 13000; s++) sph.step(dt);           // spin up well past the initial-fill splash so recording opens on settled, established flow

    gfx::App app(W, H, "waterfall", headless);
    gfx::FluidSurf fs; fs.init(W, H);
    fs.setLight(glm::vec3(-0.45f, -0.8f, -0.4f));
    gfx::Mesh box = gfx::makeBox(), plane = gfx::makePlane(40, 1);

    // foam / mist billboards for the aerated whitewater (fast, low-density particles)
    GLuint pFoam = gfx::program(R"(#version 330 core
layout(location=0) in vec3 aC; layout(location=1) in vec2 aUV; layout(location=2) in float aS;
uniform mat4 uV,uP; uniform vec3 uRight,uUp; out vec2 vUV;
void main(){ vUV=aUV; vec3 wp=aC+(aUV.x*uRight+aUV.y*uUp)*aS; gl_Position=uP*uV*vec4(wp,1.0); })",
        R"(#version 330 core
in vec2 vUV; out vec4 frag; uniform float uAlpha;
void main(){ float r=length(vUV); if(r>1.0) discard; float a=smoothstep(1.0,0.12,r)*uAlpha; frag=vec4(vec3(0.82,0.89,0.93),a); })");
    GLuint foamVAO, foamVBO; glGenVertexArrays(1, &foamVAO); glGenBuffers(1, &foamVBO);
    glBindVertexArray(foamVAO); glBindBuffer(GL_ARRAY_BUFFER, foamVBO);
    glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2); glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(5 * sizeof(float)));
    glBindVertexArray(0);

    gfx::OrbitCamera cam; cam.target = glm::vec3(0.0f, 0.95f, 0.15f); cam.dist = 5.2f; cam.yaw = 0.62f; cam.pitch = 0.15f;
    glm::mat4 proj = glm::perspective(glm::radians(43.0f), (float)W / H, 0.05f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    auto sstep = [](float a, float b, float x) { float t = (x - a) / (b - a); t = t < 0 ? 0 : (t > 1 ? 1 : t); return t * t * (3 - 2 * t); };
    std::vector<float4> host; std::vector<float> dens, surf, foam;
    const float rho0 = sph.p.rho0;
    const glm::vec3 rockA(0.31, 0.29, 0.26), rockB(0.37, 0.35, 0.31);
    const float uvq[6][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, -1}, {1, 1}, {-1, 1}};
    auto renderFrame = [&]() {
        sph.downloadPosSpeed(host); sph.downloadDens(dens);
        surf.clear(); foam.clear();
        for (size_t i = 0; i < host.size(); i++) { auto& q = host[i];
            float aer = 1.0f - sstep(0.25f * rho0, 0.75f * rho0, dens[i]);       // 0 dense bulk … 1 aerated/sparse
            float spd = sstep(2.6f, 6.0f, q.w);
            float key = std::min(1.0f, 0.55f * aer + 0.55f * spd);              // whitewater fraction (graded, not blanket white)
            if (dens[i] > 0.36f * rho0)                            // the water body — incl. the fall — as a blue-green→white surface
                { surf.push_back(q.x); surf.push_back(q.y); surf.push_back(q.z); surf.push_back(key); }
            if (aer > 0.62f && q.w > 4.2f)                         // only the fastest, most aerated water (the plunge) → light foam-mist
                for (auto& c : uvq) { float v[6] = {q.x, q.y, q.z, c[0], c[1], 0.024f}; foam.insert(foam.end(), v, v + 6); }
        }
        glm::vec3 eye = cam.eye(), fwd = glm::normalize(cam.target - eye);
        glm::vec3 R = glm::normalize(glm::cross(fwd, glm::vec3(0, 1, 0))), U = glm::cross(R, fwd);
        fs.setCamera(cam.view(), proj, eye);
        fs.beginScene(glm::vec3(0.55, 0.66, 0.80));
        fs.drawMesh(plane, glm::mat4(1), glm::vec3(0.34, 0.33, 0.31), 0.95f, 1);
        for (size_t b = 0; b < boxes.size(); b++) fs.drawMesh(box, boxModel(boxes[b]), (b & 1) ? rockB : rockA, 0.85f, 0);
        // whitewater foam-mist into the scene colour buffer (depth-tested against the rock), under the water surface
        glBindFramebuffer(GL_FRAMEBUFFER, fs.sceneFBO); glViewport(0, 0, W, H);
        GLenum c0 = GL_COLOR_ATTACHMENT0; glDrawBuffers(1, &c0);
        glUseProgram(pFoam); gfx::setM4(pFoam, "uV", cam.view()); gfx::setM4(pFoam, "uP", proj);
        gfx::setV3(pFoam, "uRight", R); gfx::setV3(pFoam, "uUp", U); gfx::setF(pFoam, "uAlpha", 0.07f);
        glEnable(GL_DEPTH_TEST); glDepthMask(GL_FALSE); glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(foamVAO); glBindBuffer(GL_ARRAY_BUFFER, foamVBO);
        glBufferData(GL_ARRAY_BUFFER, foam.size() * sizeof(float), foam.data(), GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(foam.size() / 6));
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
        GLenum bufs2[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1}; glDrawBuffers(2, bufs2);
        fs.surface(surf, 4, 0.017f);
        fs.blur(18, 0.055f);
        fs.composite(1, glm::vec3(0.05, 0.28, 0.40), glm::vec3(0.72, 0.83, 0.89));   // deep blue-green water → aerated whitewater
    };
    auto step = [&]() { for (int s = 0; s < SUB; s++) sph.step(dt); };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw += 0.30f / frames; renderFrame();     // gentle drift for parallax
            char q[512]; std::snprintf(q, sizeof(q), "%s_%04d.png", video, f); fs.screenshot(q); step(); }
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
