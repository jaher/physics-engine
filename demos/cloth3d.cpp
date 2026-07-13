// 3D cloth demo: a square sheet drapes over a sphere and folds onto the floor,
// simulated with Verlet + iterative distance constraints (phys::Cloth), gently
// disturbed by wind. Rendered as two-sided woven fabric with PBR + shadows.
//   ./cloth3d               interactive (drag to orbit, scroll to zoom)
//   ./cloth3d --shot out.png [frames]
#include "phys/phys.h"
#include "phys/cloth.h"
#include "common/gfx.h"
#include <vector>
#include <cstring>
#include <cmath>
using namespace phys;

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 320;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; frames = 260; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    const int N = 74; const real spacing = 0.055;
    real span = spacing * (N - 1);
    const Vector3 sphC(0, 1.05, 0); const real sphR = 1.05;      // sphere resting on the floor
    Cloth cloth(N, N, spacing, Vector3(-span / 2, 2.65, -span / 2), Vector3(1, 0, 0), Vector3(0, 0, 1));
    cloth.iterations = 28; cloth.damping = 0.99; cloth.friction = 0.5;
    cloth.spheres.push_back({sphC, sphR});
    cloth.groundY = 0.02;

    gfx::App app(W, H, "cloth", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);
    gfx::Mesh clothMesh;
    std::vector<unsigned> idx(cloth.tris.begin(), cloth.tris.end());
    std::vector<float> verts(cloth.pos.size() * 8);
    auto buildVerts = [&]() {
        for (size_t i = 0; i < cloth.pos.size(); i++) {
            int x = i % N, y = i / N; float* v = &verts[i * 8];
            v[0] = cloth.pos[i].x; v[1] = cloth.pos[i].y; v[2] = cloth.pos[i].z;
            v[3] = cloth.normal[i].x; v[4] = cloth.normal[i].y; v[5] = cloth.normal[i].z;
            v[6] = (float)x / (N - 1); v[7] = (float)y / (N - 1);
        }
    };
    buildVerts(); clothMesh.upload(verts, idx, true);

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 0.95, 0); cam.dist = 6.0f; cam.yaw = 0.7f; cam.pitch = 0.34f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    glm::vec3 clothColour(0.62, 0.12, 0.14), sphColour(0.28, 0.30, 0.34), groundColour(0.60, 0.58, 0.54);
    double t = 0;
    auto simulate = [&](double time) {
        // very light breeze — just enough to ripple the folds
        Vector3 wind(std::sin(time * 0.5) * 0.12, 0.0, std::sin(time * 1.1) * 0.12);
        for (int s = 0; s < 4; s++) cloth.step(1.0 / 240, wind);
        buildVerts(); clothMesh.updateVerts(verts);
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 1.0, 0), 5.0f);
        glm::mat4 sphM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(sphC.x, sphC.y, sphC.z)), glm::vec3((float)sphR));
        r.beginShadow();
        r.shadowDraw(sphere, sphM); r.shadowDraw(clothMesh, glm::mat4(1));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), groundColour, 0.95f, 0.0f);
        r.drawPBR(r.pSolid, sphere, sphM, sphColour, 0.4f, 0.05f);
        glDisable(GL_CULL_FACE);
        r.drawPBR(r.pCloth, clothMesh, glm::mat4(1), clothColour, 0.8f, 0.0f);
        r.endScene();
    };

    if (video) {                                          // render an image sequence of the fall + drape
        for (int f = 0; f < frames; f++) {
            cam.yaw = 0.55f + f * (0.55f / frames);        // slow cinematic orbit
            renderFrame();
            char path[512]; std::snprintf(path, sizeof(path), "%s_%04d.png", video, f);
            r.screenshot(path);
            simulate(t); t += 1.0 / 60;
        }
        std::printf("wrote %d frames to %s_*.png\n", frames, video); return 0;
    }
    if (headless) {
        for (int i = 0; i < frames; i++) { simulate(t); t += 1.0 / 60; }
        renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0;
    }
    double lx = 0, ly = 0; bool dragging = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && dragging) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; dragging = down;
        simulate(t); t += 1.0 / 60;
        renderFrame(); r.present(); app.poll();
    }
    return 0;
}
