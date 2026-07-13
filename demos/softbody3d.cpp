// Tetrahedral co-rotational FEM soft bodies (phys::SoftBody). Several jelly cubes
// of different stiffness drop onto the floor, squash, and wobble back — the
// per-element polar-decomposition rotation keeps the deformation rotation-stable.
//   ./softbody3d [--shot out.png [f]] | --video f [n]
#include "phys/phys.h"
#include "phys/softbody.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <array>
#include <algorithm>
using namespace phys;

// boundary triangles of a tet mesh, wound outward (faces shared by one tet only)
static std::vector<std::array<int, 3>> surfaceTris(const SoftBody& sb) {
    std::map<std::array<int, 3>, std::pair<std::array<int, 3>, int>> faces;
    const int fv[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    const int opp[4] = {3, 2, 1, 0};
    for (const auto& t : sb.tets) for (int f = 0; f < 4; f++) {
        int a = t.n[fv[f][0]], b = t.n[fv[f][1]], c = t.n[fv[f][2]], o = t.n[opp[f]];
        Vector3 n = (sb.rest[b] - sb.rest[a]) % (sb.rest[c] - sb.rest[a]);
        if (n * (sb.rest[a] - sb.rest[o]) < 0) std::swap(b, c);        // wind outward
        std::array<int, 3> key = {a, b, c}; std::sort(key.begin(), key.end());
        auto it = faces.find(key);
        if (it == faces.end()) faces[key] = {{a, b, c}, 1}; else it->second.second++;
    }
    std::vector<std::array<int, 3>> out;
    for (auto& kv : faces) if (kv.second.second == 1) out.push_back(kv.second.first);
    return out;
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int N = 5; const real sz = 0.10;                            // 5³ cells, edge 0.5

    struct Body { SoftBody sb; glm::vec3 col; };
    std::vector<Body> bodies;
    struct Spec { real x, y, youngs; glm::vec3 col; };
    Spec specs[] = {
        {-1.5, 1.6, 2.2e4, {0.86, 0.30, 0.34}}, {-0.5, 2.3, 6e4, {0.90, 0.62, 0.20}},
        {0.6, 1.5, 1.2e5, {0.30, 0.62, 0.42}},  {1.6, 2.7, 3.5e4, {0.32, 0.48, 0.78}}};
    for (auto& s : specs) {
        Body b; b.col = s.col;
        b.sb = SoftBody::makeBox(N, N, N, sz, s.youngs, 0.30, 1000, Vector3(s.x - N * sz * 0.5, s.y, -N * sz * 0.5));
        b.sb.gravity = Vector3(0, -9.81, 0); b.sb.groundY = 0.0; b.sb.damping = 0.996; b.sb.restitution = 0.1;
        bodies.push_back(std::move(b));
    }
    auto tris = surfaceTris(bodies[0].sb);                            // identical topology for all

    gfx::App app(W, H, "soft bodies", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh plane = gfx::makePlane(40, 1);
    GLuint vao, vbo; glGenVertexArrays(1, &vao); glGenBuffers(1, &vbo);
    glBindVertexArray(vao); glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glBindVertexArray(0);
    std::vector<float> vd;
    auto buildMesh = [&](const SoftBody& sb) {
        vd.clear();
        for (auto& t : tris) {
            Vector3 a = sb.nodes[t[0]].x, b = sb.nodes[t[1]].x, c = sb.nodes[t[2]].x;
            Vector3 n = (b - a) % (c - a); if (n.squareMagnitude() > 1e-12) n.normalise();
            for (Vector3 p : {a, b, c}) vd.insert(vd.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z});
        }
        glBindBuffer(GL_ARRAY_BUFFER, vbo); glBufferData(GL_ARRAY_BUFFER, vd.size() * sizeof(float), vd.data(), GL_DYNAMIC_DRAW);
        return (GLsizei)(vd.size() / 6);
    };

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 0.6, 0); cam.dist = 4.6f; cam.yaw = 0.42f; cam.pitch = 0.22f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 0.8, 0), 4.0f);
        r.beginShadow();
        for (auto& b : bodies) { GLsizei n = buildMesh(b.sb); glBindVertexArray(vao);
            gfx::setM4(r.pDepth, "uModel", glm::mat4(1)); glDrawArrays(GL_TRIANGLES, 0, n); }
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.24, 0.25, 0.28), 0.95f, 0.0f);
        glDisable(GL_CULL_FACE);
        for (auto& b : bodies) { GLsizei n = buildMesh(b.sb);
            r.useProgram(r.pSolid);
            gfx::setM4(r.pSolid, "uModel", glm::mat4(1)); gfx::setV3(r.pSolid, "uAlbedo", b.col);
            gfx::setF(r.pSolid, "uRough", 0.35f); gfx::setF(r.pSolid, "uMetal", 0.0f); gfx::setF(r.pSolid, "uAO", 1.0f);
            glBindVertexArray(vao); glDrawArrays(GL_TRIANGLES, 0, n); }
        r.endScene();
    };
    auto step = [&]() { for (int s = 0; s < 5; s++) for (auto& b : bodies) b.sb.step(1.0 / 300); };

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
