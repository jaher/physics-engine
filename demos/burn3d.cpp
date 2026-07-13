// Burning paper: a sheet hangs from a rod and is lit along its bottom edge. A
// fire-propagation cellular automaton (phys::BurningPaper) drives a glowing
// ragged burn front that climbs upward, blackening the paper to char, curling
// the hot edge, and eating a growing hole as consumed cells drop away. Smoke and
// embers rise from the flame front. ./burn3d [--shot out.png [f]] | --video f [n]
#include "phys/phys.h"
#include "phys/burn.h"
#include "common/gfx.h"
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <vector>
#include <zlib.h>
using namespace phys;

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int GW = 116, GH = 64; const real SP = 0.0115;

    BurningPaper paper;
    paper.build(GW, GH, SP, Vector3(-GW * SP * 0.5, 2.65, 0), Vector3(1, 0, 0), Vector3(0, -1, 0));
    paper.spread = 0.8f; paper.burnRate = 0.5f; paper.upBias = 1.9f; paper.cool = 0.6f; paper.curlForce = 7.5f; paper.buoyancy = 0.4f; paper.ignite = 0.5f;
    paper.flutter = 0.85f; paper.flutterSpeed = 3.2f;   // light breeze — the map undulates as it burns
    for (int x = 0; x < GW; x++) paper.pin(x, 0);                 // clipped along the top
    unsigned sd = 7; auto rnd = [&]() { sd = sd * 1103515245u + 12345u; return (float)(((sd >> 16) & 0x7fff) / 32767.0); };
    for (int y = GH - 5; y < GH; y++) for (int x = GW - 5; x < GW; x++)   // light the bottom-right corner
        paper.ignite_at(x, y, 1.1f + 0.3f * rnd());

    gfx::App app(W, H, "burning paper", headless);
    gfx::Renderer r; r.init(W, H);

    // the sheet is a scan of the 1900 Larousse world planisphere (public domain),
    // packed as zlib-compressed RGB: "PMZ1" + u32 w + u32 h + deflate(rgb)
    GLuint mapTex = 0;
    {
        const char* path = "demos/assets/worldmap.pmz";
        FILE* fp = std::fopen(path, "rb");
        if (!fp) { std::fprintf(stderr, "cannot open %s (run from the repo root)\n", path); return 1; }
        std::fseek(fp, 0, SEEK_END); long fsz = std::ftell(fp); std::fseek(fp, 0, SEEK_SET);
        std::vector<unsigned char> buf(fsz);
        if (std::fread(buf.data(), 1, fsz, fp) != (size_t)fsz) { std::fclose(fp); return 1; }
        std::fclose(fp);
        uint32_t tw = buf[4] | buf[5] << 8 | buf[6] << 16 | (uint32_t)buf[7] << 24;
        uint32_t th = buf[8] | buf[9] << 8 | buf[10] << 16 | (uint32_t)buf[11] << 24;
        uLongf rawlen = (uLongf)tw * th * 3;
        std::vector<unsigned char> rgb(rawlen);
        if (uncompress(rgb.data(), &rawlen, buf.data() + 12, (uLong)(fsz - 12)) != Z_OK) {
            std::fprintf(stderr, "world-map inflate failed\n"); return 1; }
        glGenTextures(1, &mapTex); glBindTexture(GL_TEXTURE_2D, mapTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB8, tw, th, 0, GL_RGB, GL_UNSIGNED_BYTE, rgb.data());
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
#ifdef GL_TEXTURE_MAX_ANISOTROPY_EXT
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 8.0f);
#endif
    }
    r.lightColor = glm::vec3(1.0f, 0.96f, 0.9f) * 2.7f;          // dim room so the fire glows
    r.lightDir = glm::normalize(glm::vec3(-0.4f, -0.85f, -0.5f));
    gfx::Mesh box = gfx::makeBox(), plane = gfx::makePlane(30, 1);

    // dynamic paper mesh: pos(3) normal(3) uv+char+emissive(4)
    GLuint pvao, pvbo; glGenVertexArrays(1, &pvao); glGenBuffers(1, &pvbo);
    glBindVertexArray(pvao); glBindBuffer(GL_ARRAY_BUFFER, pvbo);
    glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2); glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)(6 * sizeof(float)));
    glBindVertexArray(0);
    std::vector<float> pv; GLsizei paperVerts = 0;
    auto buildPaper = [&]() {
        pv.clear();
        auto emit = [&](int i) {
            float u = (float)(i % GW) / (GW - 1), v = (float)(i / GW) / (GH - 1);
            float ch = std::min(1.0f, paper.charAmt[i]);
            float em = glm::clamp((paper.heat[i] - 0.42f) * 2.1f, 0.0f, 1.7f);
            const Vector3& p = paper.pos[i]; const Vector3& n = paper.normal[i];
            pv.insert(pv.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z, u, v, ch, em});
        };
        for (size_t t = 0; t + 2 < paper.tris.size(); t += 3) {
            int a = paper.tris[t], b = paper.tris[t + 1], c = paper.tris[t + 2];
            if (!paper.cellAlive(a) || !paper.cellAlive(b) || !paper.cellAlive(c)) continue;
            emit(a); emit(b); emit(c);
        }
        paperVerts = (GLsizei)(pv.size() / 10);
        glBindBuffer(GL_ARRAY_BUFFER, pvbo);
        glBufferData(GL_ARRAY_BUFFER, pv.size() * sizeof(float), pv.data(), GL_DYNAMIC_DRAW);
    };

    // particles (smoke + embers), camera-facing quads: pos(3) rgba(4) uv(2)
    struct Part { Vector3 p, v; float life, max, size, r, g, b, a; };
    std::vector<Part> smoke, ember;
    GLuint qvao, qvbo; glGenVertexArrays(1, &qvao); glGenBuffers(1, &qvbo);
    glBindVertexArray(qvao); glBindBuffer(GL_ARRAY_BUFFER, qvbo);
    glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2); glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(7 * sizeof(float)));
    glBindVertexArray(0);
    std::vector<float> qv;
    auto buildQuads = [&](std::vector<Part>& ps, glm::vec3 right, glm::vec3 up) {
        qv.clear();
        for (auto& p : ps) {
            float lf = p.life / p.max;
            float alpha = p.a * (lf < 0.2f ? lf / 0.2f : (1 - lf) / 0.8f);   // fade in/out
            glm::vec3 c(p.r, p.g, p.b);
            glm::vec3 C(p.p.x, p.p.y, p.p.z);
            glm::vec3 rr = right * p.size, uu = up * p.size;
            glm::vec3 v0 = C - rr - uu, v1 = C + rr - uu, v2 = C + rr + uu, v3 = C - rr + uu;
            float quad[6][9] = {
                {v0.x,v0.y,v0.z,c.r,c.g,c.b,alpha,0,0}, {v1.x,v1.y,v1.z,c.r,c.g,c.b,alpha,1,0}, {v2.x,v2.y,v2.z,c.r,c.g,c.b,alpha,1,1},
                {v0.x,v0.y,v0.z,c.r,c.g,c.b,alpha,0,0}, {v2.x,v2.y,v2.z,c.r,c.g,c.b,alpha,1,1}, {v3.x,v3.y,v3.z,c.r,c.g,c.b,alpha,0,1}};
            for (auto& q : quad) qv.insert(qv.end(), q, q + 9);
        }
        glBindBuffer(GL_ARRAY_BUFFER, qvbo);
        glBufferData(GL_ARRAY_BUFFER, qv.size() * sizeof(float), qv.data(), GL_DYNAMIC_DRAW);
        return (GLsizei)(qv.size() / 9);
    };

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 2.32, 0); cam.dist = 2.85f; cam.yaw = 0.20f; cam.pitch = 0.03f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    glm::mat4 rod = glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, 2.68f, 0)), glm::vec3(0.70f, 0.02f, 0.02f));
    double t = 0;
    auto spawn = [&]() {
        for (int y = 0; y < GH; y++) for (int x = 0; x < GW; x++) { int i = paper.idx(x, y);
            float f = paper.fuel[i];
            if (f < 0.12f || f > 0.8f || paper.heat[i] < paper.ignite * 0.9f) continue;   // active front ring only
            if (rnd() < 0.012f) { Vector3 pp = paper.pos[i]; pp.z += (rnd() - 0.5f) * 0.02f;
                smoke.push_back({pp, Vector3((rnd() - 0.5f) * 0.3, 0.9 + rnd() * 0.6, (rnd() - 0.5f) * 0.3), 0, 1.3f + rnd() * 0.9f, 0.025f, 0.16f, 0.15f, 0.14f, 0.11f}); }
            if (rnd() < 0.009f) { Vector3 pp = paper.pos[i];
                ember.push_back({pp, Vector3((rnd() - 0.5f) * 0.55, 1.1 + rnd() * 1.0, (rnd() - 0.5f) * 0.55), 0, 0.35f + rnd() * 0.45f, 0.009f, 4.5f, 1.4f, 0.25f, 1.0f}); }
        }
    };
    auto updateParts = [&](std::vector<Part>& ps, real dt, bool grow) {
        for (auto& p : ps) { p.life += dt; p.v.y += (grow ? 0.5 : 1.5) * dt; p.v *= 0.99;
            p.p += p.v * dt; if (grow) p.size += 0.16f * dt; }
        ps.erase(std::remove_if(ps.begin(), ps.end(), [](const Part& p) { return p.life >= p.max; }), ps.end());
    };
    auto stepPhysics = [&](real dt) {
        Vector3 wind(std::sin(t * 0.7) * 0.12, 0, 0.10 + std::sin(t * 1.1) * 0.12);   // gentle gusting breeze billows the sheet
        paper.step(dt, wind);
        spawn();
        updateParts(smoke, dt, true); updateParts(ember, dt, false);
        t += dt;
    };
    auto renderFrame = [&]() {
        buildPaper();
        r.setLightForScene(glm::vec3(0, 1.6, 0), 3.5f);
        r.beginShadow();
        r.shadowDraw(box, rod);
        gfx::setM4(r.pDepth, "uModel", glm::mat4(1)); glBindVertexArray(pvao); glDrawArrays(GL_TRIANGLES, 0, paperVerts);
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.14, 0.13, 0.13), 0.9f, 0.0f);
        r.drawPBR(r.pSolid, box, rod, glm::vec3(0.2, 0.16, 0.1), 0.7f, 0.1f);
        // paper
        glDisable(GL_CULL_FACE);
        r.useProgram(r.pPaper);
        glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, mapTex); gfx::setI(r.pPaper, "uMap", 1);
        glActiveTexture(GL_TEXTURE0);
        glBindVertexArray(pvao); glDrawArrays(GL_TRIANGLES, 0, paperVerts);
        // particles
        glm::vec3 eye = cam.eye(), fwd = glm::normalize(cam.target - eye);
        glm::vec3 rt = glm::normalize(glm::cross(fwd, glm::vec3(0, 1, 0))), upv = glm::cross(rt, fwd);
        glUseProgram(r.pPart); gfx::setM4(r.pPart, "uView", cam.view()); gfx::setM4(r.pPart, "uProj", proj);
        glEnable(GL_BLEND); glDepthMask(GL_FALSE);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        { GLsizei n = buildQuads(smoke, rt, upv); glBindVertexArray(qvao); glDrawArrays(GL_TRIANGLES, 0, n); }
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);           // additive embers
        { GLsizei n = buildQuads(ember, rt, upv); glBindVertexArray(qvao); glDrawArrays(GL_TRIANGLES, 0, n); }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 4; s++) stepPhysics(1.0 / 240); }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames * 4; i++) stepPhysics(1.0 / 240); renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        for (int s = 0; s < 4; s++) stepPhysics(1.0 / 240);
        renderFrame(); r.present(); app.poll();
    }
    return 0;
}
