// Reduced-coordinate articulated bodies (Featherstone ABA). A chain of rigid
// links joined by revolute joints is simulated in *joint* coordinates — released
// from horizontal it whips and swings as a chaotic multi-pendulum, driven by the
// O(n) articulated-body algorithm in phys::Articulation.
//   ./articulation3d [--shot out.png [f]] | --video f [n]
#include "phys/phys.h"
#include "phys/articulation.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
#include <vector>
using namespace phys;

static glm::mat4 glmFromPhys(const Matrix4& t) {
    return glm::mat4(t.data[0], t.data[4], t.data[8], 0, t.data[1], t.data[5], t.data[9], 0,
                     t.data[2], t.data[6], t.data[10], 0, t.data[3], t.data[7], t.data[11], 1);
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int NLINK = 6; const real L = 0.62, m = 1.0;
    const Vector3 pivot(0, 4.4, 0);

    Articulation arm; arm.setGravity(Vector3(0, -9.81, 0));
    Matrix3 I; I.setInertiaTensorCoeffs(m * L * L / 12, m * L * L / 48, m * L * L / 12);   // thin rod about COM
    for (int i = 0; i < NLINK; i++) {
        Vector3 jointPos = i == 0 ? pivot : Vector3(0, -L, 0);        // world for root, parent-frame otherwise
        arm.addLink(i - 1, JointType::REVOLUTE, Vector3(0, 0, 1), m, I, Vector3(0, -L / 2, 0), jointPos);
    }
    arm.setJointState(0, 1.55, 0);                                    // release the root near-horizontal

    gfx::App app(W, H, "articulated chain", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 2.6, 0); cam.dist = 6.4f; cam.yaw = 0.16f; cam.pitch = 0.05f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    glm::mat4 bracket = glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, pivot.y + 0.12f, 0)), glm::vec3(0.5f, 0.12f, 0.14f));
    glm::mat4 postL = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-0.42f, pivot.y + 0.6f, 0)), glm::vec3(0.06f, 0.6f, 0.08f));

    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 2.6, 0), 5.5f);
        r.beginShadow();
        r.shadowDraw(box, bracket); r.shadowDraw(box, postL);
        for (int i = 0; i < NLINK; i++) {
            glm::mat4 f = glmFromPhys(arm.linkWorld(i));
            r.shadowDraw(box, glm::scale(glm::translate(f, glm::vec3(0, -L / 2, 0)), glm::vec3(0.045f, L / 2, 0.045f)));
        }
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.20, 0.21, 0.24), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, bracket, glm::vec3(0.16, 0.16, 0.18), 0.6f, 0.3f);
        r.drawPBR(r.pSolid, box, postL, glm::vec3(0.16, 0.16, 0.18), 0.6f, 0.3f);
        for (int i = 0; i < NLINK; i++) {
            glm::mat4 f = glmFromPhys(arm.linkWorld(i));
            float t = (float)i / (NLINK - 1);
            glm::vec3 col = glm::mix(glm::vec3(0.86, 0.36, 0.16), glm::vec3(0.18, 0.42, 0.72), t);   // warm→cool along the chain
            r.drawPBR(r.pSolid, box, glm::scale(glm::translate(f, glm::vec3(0, -L / 2, 0)), glm::vec3(0.045f, L / 2, 0.045f)), col, 0.35f, 0.2f);
            r.drawPBR(r.pSolid, sphere, glm::scale(f, glm::vec3(0.075f)), glm::vec3(0.7, 0.72, 0.75), 0.25f, 0.85f);   // joint
        }
        r.endScene();
    };
    auto step = [&]() { for (int s = 0; s < 4; s++) arm.step(1.0 / 240); };

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
