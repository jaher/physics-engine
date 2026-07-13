// 3D hair demo: ~2600 strands grown from a scalp, each a Verlet chain with
// segment + bend constraints colliding with the head (phys::Hair), swaying in a
// breeze. Rendered as camera-facing ribbons with Kajiya-Kay anisotropic shading
// (two shifted specular lobes) for that characteristic strand sheen.
//   ./hair3d               interactive (drag to orbit, scroll to zoom)
//   ./hair3d --shot out.png [frames]
#include "phys/phys.h"
#include "phys/hair.h"
#include "common/gfx.h"
#include <vector>
#include <cstring>
#include <cmath>
using namespace phys;

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 160;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; frames = 300; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    const Vector3 headC(0, 1.75, 0); const real headR = 0.95;
    Hair hair; hair.headCentre = headC; hair.headRadius = headR;
    hair.gravity = Vector3(0, -9.81, 0); hair.damping = 0.965; hair.iterations = 16; hair.bendStiffness = 0.72;
    hair.turbulence = 2.6;                                // per-strand fluttering in the wind

    const int SEG = 15;                                   // points per strand = SEG (14 segments)
    // scalp roots via a Fibonacci sphere, keeping the top/back/sides (no face)
    unsigned seed = 12345;
    auto rnd = [&]() { seed = seed * 1103515245u + 12345u; return ((seed >> 16) & 0x7fff) / 32767.0; };
    int candidates = 8200;
    const double PHI = M_PI * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < candidates; i++) {
        double y = 1.0 - (i / (double)(candidates - 1)) * 2.0;
        double rad = std::sqrt(std::max(0.0, 1.0 - y * y));
        double phi = i * PHI;
        double x = std::cos(phi) * rad, z = std::sin(phi) * rad;
        Vector3 nrm(x, y, z);
        if (nrm.y < -0.05) continue;                                   // below the crown
        if (nrm.z > 0.42 && nrm.y < 0.55) continue;                    // leave the face bare
        Vector3 root = headC + nrm * (headR * 1.001);
        real len = 1.25 + rnd() * 0.5;
        real seg = len / (SEG - 1);
        real shade = (real)(rnd() * 2.0 - 1.0);
        hair.addStrand(root, nrm, SEG - 1, seg, 0.008, 0.0028, shade);
    }
    int S = (int)hair.strands.size();

    // Clumping: render K child strands around each simulated guide (offset in the
    // root's tangent plane, fanning out toward the tip). Dense hair at guide cost.
    const int K = 7;
    struct Child { Vector3 off; float shade; };
    std::vector<Child> children(S * K);
    { unsigned cs = 9871;
      auto cr = [&]() { cs = cs * 1103515245u + 12345u; return (double)((cs >> 16) & 0x7fff) / 32767.0; };
      for (int s = 0; s < S; s++) {
          Vector3 nrm = hair.strands[s].pos[0] - headC; nrm.normalise();
          Vector3 up(0, 1, 0); if (real_abs(nrm * up) > 0.95) up = Vector3(1, 0, 0);
          Vector3 u = nrm % up; u.normalise(); Vector3 v = nrm % u;
          for (int c = 0; c < K; c++) {
              double a = cr() * 6.2831853, mag = 0.004 + cr() * 0.028;
              Vector3 off = (u * std::cos(a) + v * std::sin(a)) * mag;
              children[s * K + c] = {off, (float)((cr() * 2 - 1) * 0.22)};
          }
      }
    }
    const int R = S * K;                                  // total rendered strands
    std::printf("%d guide strands, %d rendered strands\n", S, R);

    gfx::App app(W, H, "hair", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);

    // custom ribbon buffer: pos(3) tangent(3) tv(3)=(t,side,shade)
    GLuint vao, vbo, ebo; glGenVertexArrays(1, &vao); glGenBuffers(1, &vbo); glGenBuffers(1, &ebo);
    std::vector<unsigned> idx; idx.reserve((size_t)R * (SEG - 1) * 6);
    for (int rs = 0; rs < R; rs++) { unsigned base = rs * SEG * 2;
        for (int i = 0; i < SEG - 1; i++) { unsigned a = base + 2 * i, b = a + 1, c = a + 2, d = a + 3;
            idx.insert(idx.end(), {a, c, b, b, c, d}); } }
    glBindVertexArray(vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(unsigned), idx.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2); glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glBindVertexArray(0);
    std::vector<float> verts((size_t)R * SEG * 2 * 9);

    auto buildRibbons = [&](const Vector3& camPos) {
        for (int s = 0; s < S; s++) { auto& st = hair.strands[s]; int P = (int)st.pos.size();
            for (int c = 0; c < K; c++) { int rs = s * K + c; const Child& ch = children[rs];
                for (int i = 0; i < P; i++) {
                    real t = (real)i / (P - 1);
                    int im = i > 0 ? i - 1 : 0, ip = i < P - 1 ? i + 1 : P - 1;
                    Vector3 p  = st.pos[i]  + ch.off * t;
                    Vector3 pa = st.pos[im] + ch.off * ((real)im / (P - 1));
                    Vector3 pb = st.pos[ip] + ch.off * ((real)ip / (P - 1));
                    Vector3 tan = pb - pa; if (tan.squareMagnitude() < 1e-9) tan = Vector3(0, -1, 0); tan.normalise();
                    Vector3 toCam = camPos - p; toCam.normalise();
                    Vector3 side = tan % toCam; if (side.squareMagnitude() < 1e-9) side = Vector3(1, 0, 0); side.normalise();
                    real hw = (st.widthRoot * (1 - t) + st.widthTip * t) * (real)0.5;
                    Vector3 v0 = p - side * hw, v1 = p + side * hw;
                    float sh = (float)st.shade + ch.shade; if (sh > 1) sh = 1; if (sh < -1) sh = -1;
                    int b = (rs * SEG + i) * 2 * 9;
                    float rec0[9] = {(float)v0.x, (float)v0.y, (float)v0.z, (float)tan.x, (float)tan.y, (float)tan.z, (float)t, -1.f, sh};
                    float rec1[9] = {(float)v1.x, (float)v1.y, (float)v1.z, (float)tan.x, (float)tan.y, (float)tan.z, (float)t, 1.f, sh};
                    for (int k = 0; k < 9; k++) { verts[b + k] = rec0[k]; verts[b + 9 + k] = rec1[k]; }
                }
            }
        }
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(), GL_DYNAMIC_DRAW);
    };

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.35, 0); cam.dist = 4.6f; cam.yaw = 0.35f; cam.pitch = 0.12f;
    glm::mat4 proj = glm::perspective(glm::radians(40.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    glm::vec3 headColour(0.72, 0.56, 0.47), rootCol(0.07, 0.04, 0.022), tipCol(0.32, 0.18, 0.085);
    double t = 0;
    glm::mat4 headM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(headC.x, headC.y, headC.z)), glm::vec3((float)headR));

    auto simulate = [&](double time) {
        // a gusting wind blowing from the left, streaming the hair to the right
        double gust = 5.0 + 3.2 * std::sin(time * 0.7) + 1.8 * std::sin(time * 1.9 + 1.0);
        Vector3 dir(1.0, 0.18, 0.35); dir.normalise();
        Vector3 wind = dir * gust;
        for (int s = 0; s < 3; s++) hair.step(1.0 / 180, wind);
    };
    auto renderFrame = [&]() {
        glm::vec3 camPos = cam.eye();
        buildRibbons(Vector3(camPos.x, camPos.y, camPos.z));
        r.setLightForScene(glm::vec3(0, 1.4, 0), 3.6f);
        r.beginShadow(); r.shadowDraw(sphere, headM); r.endShadow();
        r.beginScene(cam.view(), proj, camPos);
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.58, 0.56, 0.53), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, sphere, headM, headColour, 0.6f, 0.0f);
        // hair
        r.useProgram(r.pHair);
        gfx::setV3(r.pHair, "uRootColor", rootCol); gfx::setV3(r.pHair, "uTipColor", tipCol);
        glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(vao); glDrawElements(GL_TRIANGLES, (GLsizei)idx.size(), GL_UNSIGNED_INT, 0);
        glDisable(GL_BLEND);
        r.endScene();
    };

    if (video) {                                          // hair blowing in the wind
        hair.turbulence = 0; for (int i = 0; i < 150; i++) hair.step(1.0 / 180, Vector3());   // settle calmly first
        hair.turbulence = 2.6;
        for (int f = 0; f < frames; f++) {
            cam.yaw = 0.35f + f * (0.35f / frames);        // gentle drift
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
