// Spring harmonics: six coil springs strung between two posts, each excited in a
// pure normal mode (1…6 antinodes) of a coupled spring-mass chain (phys::SpringChain).
// They ring as standing waves — the harmonic series — the higher modes visibly
// faster, exactly as a plucked string's overtones. Each coil is a helix swept as a
// metal tube that follows the live mass positions.
//   ./harmonics3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "phys/springs.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
#include <memory>
#include <vector>
using namespace phys;

// A coil spring: a helix of radius coilR wound turnsPerUnit times along the chain
// centreline, swept as a round wire tube of radius wireR. Rebuilt every frame from
// the live node positions so the standing wave shows in the coils.
struct CoilMesh {
    gfx::Mesh mesh; std::vector<float> verts; std::vector<unsigned> idx; bool init = false;
    static Vector3 catmull(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, real t) {
        real t2 = t * t, t3 = t2 * t;
        return (p1 * 2 + (p2 - p0) * t + (p0 * 2 - p1 * 5 + p2 * 4 - p3) * t2 + (p1 * 3 - p0 - p2 * 3 + p3) * t3) * (real)0.5;
    }
    void build(const std::vector<Vector3>& nodes, real coilR, real wireR, real turnsPerUnit) {
        const int SUB = 8, RING = 8;
        int n = (int)nodes.size();
        std::vector<Vector3> c;                                   // smoothed centreline
        for (int i = 0; i + 1 < n; i++) {
            const Vector3& p0 = nodes[i > 0 ? i - 1 : 0]; const Vector3& p1 = nodes[i];
            const Vector3& p2 = nodes[i + 1]; const Vector3& p3 = nodes[i + 2 < n ? i + 2 : n - 1];
            for (int s = 0; s < SUB; s++) c.push_back(catmull(p0, p1, p2, p3, (real)s / SUB));
        }
        c.push_back(nodes.back());
        int C = (int)c.size();

        // wind a helix around the centreline (parallel-transported frame)
        std::vector<Vector3> hel(C);
        {
            Vector3 tan0 = (c[1] - c[0]).unit();
            Vector3 up = real_abs(tan0.y) < 0.95 ? Vector3(0, 1, 0) : Vector3(1, 0, 0);
            Vector3 nrm = (tan0 % up).unit(); real arc = 0;
            for (int i = 0; i < C; i++) {
                Vector3 tan = (c[i < C - 1 ? i + 1 : i] - c[i > 0 ? i - 1 : i]).unit();
                nrm = (nrm - tan * (nrm * tan)).unit(); Vector3 bin = (tan % nrm).unit();
                if (i > 0) arc += (c[i] - c[i - 1]).magnitude();
                real phi = 2 * real_pi * turnsPerUnit * arc;
                hel[i] = c[i] + (nrm * real_cos(phi) + bin * real_sin(phi)) * coilR;
            }
        }
        // sweep a round tube along the helix
        int R = C; verts.assign((size_t)R * RING * 6, 0.0f);
        std::vector<Vector3> grid((size_t)R * RING);
        Vector3 t0 = (hel[1] - hel[0]).unit();
        Vector3 up = real_abs(t0.y) < 0.95 ? Vector3(0, 1, 0) : Vector3(1, 0, 0);
        Vector3 nrm = (t0 % up).unit();
        for (int i = 0; i < R; i++) {
            Vector3 tan = (hel[i < R - 1 ? i + 1 : i] - hel[i > 0 ? i - 1 : i]).unit();
            nrm = (nrm - tan * (nrm * tan)).unit(); Vector3 bin = (tan % nrm).unit();
            for (int k = 0; k < RING; k++) {
                real th = 2 * real_pi * k / RING;
                grid[(size_t)i * RING + k] = hel[i] + (nrm * real_cos(th) + bin * real_sin(th)) * wireR;
            }
        }
        idx.clear();
        for (int i = 0; i + 1 < R; i++) for (int k = 0; k < RING; k++) {
            unsigned a = i * RING + k, b = i * RING + (k + 1) % RING, d = (i + 1) * RING + k, e = (i + 1) * RING + (k + 1) % RING;
            idx.insert(idx.end(), {a, d, b, b, d, e});
        }
        std::vector<Vector3> vn((size_t)R * RING);
        for (size_t f = 0; f < idx.size(); f += 3) {
            Vector3 fn = (grid[idx[f + 1]] - grid[idx[f]]) % (grid[idx[f + 2]] - grid[idx[f]]);
            vn[idx[f]] += fn; vn[idx[f + 1]] += fn; vn[idx[f + 2]] += fn;
        }
        for (size_t v = 0; v < grid.size(); v++) { Vector3 nn = vn[v]; if (nn.squareMagnitude() > 1e-12) nn.normalise();
            float* o = &verts[v * 6];
            o[0] = grid[v].x; o[1] = grid[v].y; o[2] = grid[v].z; o[3] = nn.x; o[4] = nn.y; o[5] = nn.z; }
        if (!init) { mesh.upload(verts, idx); init = true; } else mesh.updateVerts(verts);
    }
};

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int NMODES = 6, P = 30;
    const real L = 5.0, spacing = L / (P + 1);
    const real k = 22.0, m = 0.02, rest = spacing * 0.45;         // pre-tensioned → crisp transverse modes
    const real amp = 0.25, laneGap = 0.80, yTop = 4.5, x0 = -L / 2;

    std::vector<std::unique_ptr<SpringChain>> chains;
    for (int n = 1; n <= NMODES; n++) {
        auto c = std::make_unique<SpringChain>();
        real y = yTop - (n - 1) * laneGap;
        c->build(P, spacing, k, m, Vector3(x0, y, 0), Vector3(1, 0, 0), rest, /*damping*/0.9997);
        c->setMode(n, amp, Vector3(0, 1, 0));                     // transverse standing wave, n antinodes
        chains.push_back(std::move(c));
    }

    gfx::App app(W, H, "spring harmonics", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), plane = gfx::makePlane(30, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 2.5, 0); cam.dist = 8.2f; cam.yaw = 0.30f; cam.pitch = 0.13f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    // metallic tint per harmonic, fundamental → 6th
    glm::vec3 tint[NMODES] = {
        {0.80f, 0.48f, 0.22f}, {0.78f, 0.60f, 0.28f}, {0.85f, 0.70f, 0.24f},
        {0.60f, 0.62f, 0.66f}, {0.70f, 0.73f, 0.78f}, {0.55f, 0.68f, 0.72f}};

    // two posts the springs are strung between
    real yBot = yTop - (NMODES - 1) * laneGap;
    real postH = (yTop - yBot) * 0.5f + 0.35f, postY = (yTop + yBot) * 0.5f;
    glm::mat4 postL = glm::scale(glm::translate(glm::mat4(1), glm::vec3(x0 - 0.06f, postY, 0)), glm::vec3(0.05f, postH, 0.05f));
    glm::mat4 postR = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-x0 + 0.06f, postY, 0)), glm::vec3(0.05f, postH, 0.05f));

    std::vector<CoilMesh> coils(NMODES);
    std::vector<Vector3> nodesBuf;
    auto rebuild = [&](int n) {
        nodesBuf.clear();
        for (auto& p : chains[n]->node) nodesBuf.push_back(p.getPosition());
        coils[n].build(nodesBuf, 0.05, 0.013, 4.0);
    };

    auto renderFrame = [&]() {
        for (int n = 0; n < NMODES; n++) rebuild(n);
        r.setLightForScene(glm::vec3(0, 1.9, 0), 6.5f);
        r.beginShadow();
        r.shadowDraw(box, postL); r.shadowDraw(box, postR);
        for (int n = 0; n < NMODES; n++) r.shadowDraw(coils[n].mesh, glm::mat4(1));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.20, 0.21, 0.23), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, postL, glm::vec3(0.12, 0.12, 0.13), 0.5f, 0.6f);
        r.drawPBR(r.pSolid, box, postR, glm::vec3(0.12, 0.12, 0.13), 0.5f, 0.6f);
        for (int n = 0; n < NMODES; n++)
            r.drawPBR(r.pSolid, coils[n].mesh, glm::mat4(1), tint[n], 0.28f, 1.0f);
        r.endScene();
    };
    auto stepAll = [&](real dt) { for (auto& c : chains) c->step(dt); };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 0.30f + 0.35f * std::sin(f * 6.2831853 / frames); renderFrame();
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
