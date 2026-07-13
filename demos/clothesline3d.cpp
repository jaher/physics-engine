// Clothing hanging on a line: several cloth sheets pinned along a sagging string
// strung between two posts, hanging under gravity and swaying in a breeze
// (phys::Cloth). Two-sided woven PBR + shadows.
//   ./clothesline3d               interactive
//   ./clothesline3d --shot out.png [frames]   |   --video frames/f [n]
#include "phys/phys.h"
#include "phys/cloth.h"
#include "common/gfx.h"
#include <vector>
#include <cstring>
#include <cmath>
using namespace phys;

static const double TOPY = 3.1, SAG = 0.30, SPAN = 3.35;
static double catY(double x) { return TOPY - SAG * (1.0 - (x / SPAN) * (x / SPAN)); }   // sagging line

struct Garment {
    Cloth cloth; int N, M; glm::vec3 colour; gfx::Mesh mesh; std::vector<float> verts;
    Garment(double cx, double halfW, int N, int M, double sp, double z, glm::vec3 col)
        : cloth(N, M, sp, Vector3(cx - halfW, catY(cx - halfW), z), Vector3(1, 0, 0), Vector3(0, -1, 0)),
          N(N), M(M), colour(col) {
        cloth.iterations = 24; cloth.damping = 0.99; cloth.groundY = 0.02; cloth.friction = 0.4;
        for (int x = 0; x < N; x++) {                                   // pin the top edge to the string
            double px = cx - halfW + sp * x;
            cloth.setPin(x, 0, Vector3(px, catY(px), z));
        }
        verts.resize((size_t)N * M * 8);
        std::vector<unsigned> idx(cloth.tris.begin(), cloth.tris.end());
        build(); mesh.upload(verts, idx, true);
    }
    void build() {
        for (int i = 0; i < N * M; i++) { int x = i % N, y = i / N; float* v = &verts[i * 8];
            v[0] = cloth.pos[i].x; v[1] = cloth.pos[i].y; v[2] = cloth.pos[i].z;
            v[3] = cloth.normal[i].x; v[4] = cloth.normal[i].y; v[5] = cloth.normal[i].z;
            v[6] = (float)x / (N - 1) * 3; v[7] = (float)y / (M - 1) * 3; }
    }
    void step(real dt, const Vector3& wind) { for (int s = 0; s < 4; s++) cloth.step(dt, wind); build(); mesh.updateVerts(verts); }
};

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    gfx::App app(W, H, "clothesline", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh sphere = gfx::makeSphere(), box = gfx::makeBox(), plane = gfx::makePlane(40, 1);

    std::vector<Garment*> garments;
    garments.push_back(new Garment(-2.15, 0.95, 40, 46, 0.05, 0.0, glm::vec3(0.85, 0.85, 0.88)));   // white sheet
    garments.push_back(new Garment(0.0,  0.85, 36, 52, 0.05, 0.0, glm::vec3(0.18, 0.34, 0.62)));    // blue
    garments.push_back(new Garment(2.05, 0.90, 38, 42, 0.05, 0.0, glm::vec3(0.66, 0.20, 0.18)));    // red

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.6, 0); cam.dist = 8.0f; cam.yaw = 0.15f; cam.pitch = 0.10f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    // string rendered as a row of small beads along the catenary
    std::vector<glm::mat4> beads;
    for (int i = 0; i <= 90; i++) { double x = -SPAN + 2 * SPAN * i / 90.0;
        beads.push_back(glm::scale(glm::translate(glm::mat4(1), glm::vec3(x, catY(x), 0)), glm::vec3(0.02f))); }
    glm::mat4 poleL = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-SPAN, catY(-SPAN) * 0.5f, 0)), glm::vec3(0.06f, catY(-SPAN) * 0.5f, 0.06f));
    glm::mat4 poleR = glm::scale(glm::translate(glm::mat4(1), glm::vec3(SPAN, catY(SPAN) * 0.5f, 0)), glm::vec3(0.06f, catY(SPAN) * 0.5f, 0.06f));

    double t = 0;
    auto simulate = [&](double time) {
        Vector3 wind(std::sin(time * 0.7) * 0.45, 0.0, 0.7 + std::sin(time * 1.5) * 0.5);
        for (auto* g : garments) g->step(1.0 / 240, wind);
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 1.6, 0), 6.0f);
        r.beginShadow();
        r.shadowDraw(box, poleL); r.shadowDraw(box, poleR);
        for (auto* g : garments) r.shadowDraw(g->mesh, glm::mat4(1));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.55, 0.56, 0.52), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, poleL, glm::vec3(0.35, 0.24, 0.14), 0.8f, 0.0f);
        r.drawPBR(r.pSolid, box, poleR, glm::vec3(0.35, 0.24, 0.14), 0.8f, 0.0f);
        for (auto& m : beads) r.drawPBR(r.pSolid, sphere, m, glm::vec3(0.1, 0.09, 0.08), 0.7f, 0.0f);
        glDisable(GL_CULL_FACE);
        for (auto* g : garments) r.drawPBR(r.pCloth, g->mesh, glm::mat4(1), g->colour, 0.82f, 0.0f);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 0.05f + f * (0.32f / frames);
            renderFrame(); char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            simulate(t); t += 1.0 / 60; }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames; i++) { simulate(t); t += 1.0 / 60; } renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        simulate(t); t += 1.0 / 60; renderFrame(); r.present(); app.poll();
    }
    return 0;
}
