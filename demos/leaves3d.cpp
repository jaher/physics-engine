// Dozens of oak leaves spiralling down and gradually burning. Each leaf is a thin
// plate with tumbling-plate aerodynamics (phys::FallingLeaf) — anisotropic drag,
// reactive lift and an underdamped broadside-seeking torque give the fluttering,
// tumbling, drifting descent of a real falling leaf. Partway down each one ignites
// (staggered) and a per-leaf combustion CA eats a glowing ragged front across it,
// charring it black, curling and holing it, shedding smoke and embers — the same
// fire model as the burning paper. Consumed/landed leaves respawn at the top.
//   ./leaves3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "phys/leaf.h"
#include "common/gfx.h"
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <vector>
#include <memory>
#include <zlib.h>
using namespace phys;

// textured, two-sided leaf shader: real oak-leaf photo (texture array) + per-cell
// char (blackening) and emissive (glowing burn front); alpha cuts the silhouette
static const char* VS_LEAF = R"(#version 330 core
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aNormal; layout(location=2) in vec2 aUV; layout(location=3) in vec3 aData;
uniform mat4 uView,uProj; out vec3 vN; out vec2 vUV; out float vLayer; out float vChar; out float vEmis;
void main(){ vN=aNormal; vUV=aUV; vLayer=aData.x; vChar=aData.y; vEmis=aData.z; gl_Position=uProj*uView*vec4(aPos,1.0); })";
static const char* FS_LEAF = R"(#version 330 core
in vec3 vN; in vec2 vUV; in float vLayer; in float vChar; in float vEmis; out vec4 frag;
uniform vec3 uLightDir; uniform sampler2DArray uLeaf;
void main(){ vec4 tx=texture(uLeaf, vec3(vUV, vLayer)); if(tx.a<0.4) discard;
  vec3 N=normalize(vN); if(!gl_FrontFacing) N=-N;
  float dif=max(dot(N,-normalize(uLightDir)),0.0);
  vec3 alb=mix(tx.rgb, vec3(0.05,0.04,0.035), vChar);           // char blackens the leaf
  vec3 base=alb*(0.35+0.75*dif);
  base += vec3(3.2,1.0,0.2)*vEmis;                              // glowing burn front
  frag=vec4(base,1.0); })";

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int NLEAF = 60;

    unsigned sd = 20260709u; auto rnd = [&]() { sd = sd * 1103515245u + 12345u; return (float)(((sd >> 16) & 0x7fff) / 32767.0f); };
    auto randUnit = [&]() { Vector3 v(rnd() - 0.5f, rnd() - 0.5f, rnd() - 0.5f); real m = v.magnitude(); return m > 1e-4 ? v * (1.0 / m) : Vector3(0, 1, 0); };

    // load the real oak-leaf textures: "LEFA" + u32 count,w,h + zlib(RGBA * count layers)
    int texCount = 0, texW = 0, texH = 0; std::vector<unsigned char> texRGBA;   // all layers, RGBA
    std::vector<std::vector<unsigned char>> leafAlpha;                          // per-layer alpha (for masks)
    {
        FILE* fp = std::fopen("demos/assets/leaves.lefa", "rb");
        if (!fp) { std::fprintf(stderr, "cannot open demos/assets/leaves.lefa (run from repo root)\n"); return 1; }
        std::fseek(fp, 0, SEEK_END); long fsz = std::ftell(fp); std::fseek(fp, 0, SEEK_SET);
        std::vector<unsigned char> buf(fsz); if (std::fread(buf.data(), 1, fsz, fp) != (size_t)fsz) { std::fclose(fp); return 1; } std::fclose(fp);
        texCount = buf[4] | buf[5] << 8 | buf[6] << 16 | (int)buf[7] << 24;
        texW = buf[8] | buf[9] << 8 | buf[10] << 16 | (int)buf[11] << 24;
        texH = buf[12] | buf[13] << 8 | buf[14] << 16 | (int)buf[15] << 24;
        uLongf raw = (uLongf)texCount * texW * texH * 4; texRGBA.resize(raw);
        if (uncompress(texRGBA.data(), &raw, buf.data() + 16, (uLong)(fsz - 16)) != Z_OK) { std::fprintf(stderr, "leaf inflate failed\n"); return 1; }
        leafAlpha.resize(texCount);
        for (int l = 0; l < texCount; l++) { leafAlpha[l].resize((size_t)texW * texH);
            for (int p = 0; p < texW * texH; p++) leafAlpha[l][p] = texRGBA[((size_t)l * texW * texH + p) * 4 + 3]; }
    }

    auto spawn = [&](FallingLeaf& L, bool first) {
        L.texLayer = (int)(rnd() * texCount) % texCount;
        L.buildTextured(34, 34, 0.40, 0.40, leafAlpha[L.texLayer].data(), texW, texH, sd + 1);   // shape from the leaf's alpha
        L.pos = Vector3((rnd() - 0.5f) * 3.8, first ? (1.8 + rnd() * 5.4) : (5.9 + rnd() * 1.4), (rnd() - 0.5f) * 3.2);
        L.vel = Vector3((rnd() - 0.5f) * 0.4, -0.2f * rnd(), (rnd() - 0.5f) * 0.4);
        L.e2 = randUnit(); L.e0 = (L.e2 % Vector3(0, 1, 0.3)); if (L.e0.magnitude() < 1e-3) L.e0 = Vector3(1, 0, 0); L.e0.normalise(); L.e1 = (L.e2 % L.e0);
        L.igniteAt = 1e9; L.igniteY = 2.4 + rnd() * 2.4; L.time = 0; L.ignited = false;   // catches fire partway down
        L.burnRate = 0.17; L.spread = 0.7; L.curlAmt = 0.0; L.cool = 0.7;   // flat leaves (no out-of-plane curl)
    };
    std::vector<FallingLeaf> leaves(NLEAF);
    for (auto& L : leaves) { spawn(L, true);                               // pre-roll so they start spread through the fall
        int pre = (int)(rnd() * 1500); for (int s = 0; s < pre && !L.consumed(); s++) L.step(1.0 / 240, Vector3()); }

    gfx::App app(W, H, "burning leaves", headless);
    gfx::Renderer r; r.init(W, H);
    r.lightColor = glm::vec3(1.0f, 0.93f, 0.82f) * 3.0f;
    r.lightDir = glm::normalize(glm::vec3(-0.4f, -0.85f, -0.4f));
    gfx::Mesh plane = gfx::makePlane(60, 1);
    GLuint pLeaf = gfx::program(VS_LEAF, FS_LEAF);

    // upload the leaf textures as a GL texture array
    GLuint leafTex; glGenTextures(1, &leafTex); glBindTexture(GL_TEXTURE_2D_ARRAY, leafTex);
    glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_SRGB8_ALPHA8, texW, texH, texCount, 0, GL_RGBA, GL_UNSIGNED_BYTE, texRGBA.data());
    glGenerateMipmap(GL_TEXTURE_2D_ARRAY);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // leaf mesh VBO (pos3 normal3 uv2 data3=(layer,char,emis))
    GLuint lvao, lvbo; glGenVertexArrays(1, &lvao); glGenBuffers(1, &lvbo);
    glBindVertexArray(lvao); glBindBuffer(GL_ARRAY_BUFFER, lvbo);
    glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2); glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(3); glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(8 * sizeof(float)));
    glBindVertexArray(0);
    std::vector<float> lv; GLsizei leafVerts = 0;
    auto buildLeaves = [&]() {
        lv.clear();
        for (auto& L : leaves) {
            for (int j = 0; j + 1 < L.GH; j++) for (int i = 0; i + 1 < L.GW; i++) {
                if (!L.alive(i, j) || !L.alive(i + 1, j) || !L.alive(i, j + 1) || !L.alive(i + 1, j + 1)) continue;
                Vector3 A = L.node(i, j), B = L.node(i + 1, j), C = L.node(i, j + 1), D = L.node(i + 1, j + 1);
                Vector3 fn = ((B - A) % (C - A)); if (fn.squareMagnitude() > 1e-12) fn.normalise();
                int q[4][2] = {{i, j}, {i + 1, j}, {i, j + 1}, {i + 1, j + 1}}; Vector3 wp[4] = {A, B, C, D};
                auto vert = [&](int k) {
                    int gi = q[k][0], gj = q[k][1], c = L.idx(gi, gj);
                    float uu = (float)gi / (L.GW - 1), vv = 1.0f - (float)gj / (L.GH - 1);
                    float ch = std::min(1.0f, L.charAmt[c]);
                    float em = glm::clamp((L.heat[c] - 0.62f) * 1.8f, 0.0f, 0.11f);   // thin front (only the hottest cells glow), low amplitude
                    em *= 0.72f + 0.28f * std::sin((float)L.time * 11.0f + gi * 1.9f + gj * 2.7f);   // flicker / undulate
                    lv.insert(lv.end(), {(float)wp[k].x, (float)wp[k].y, (float)wp[k].z, (float)fn.x, (float)fn.y, (float)fn.z, uu, vv, (float)L.texLayer, ch, em});
                };
                vert(0); vert(2); vert(1); vert(1); vert(2); vert(3);
            }
        }
        leafVerts = (GLsizei)(lv.size() / 11);
        glBindBuffer(GL_ARRAY_BUFFER, lvbo); glBufferData(GL_ARRAY_BUFFER, lv.size() * sizeof(float), lv.data(), GL_DYNAMIC_DRAW);
    };

    // smoke + ember particles
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
        for (auto& p : ps) { float lf = p.life / p.max; float a = p.a * (lf < 0.2f ? lf / 0.2f : (1 - lf) / 0.8f);
            glm::vec3 c(p.r, p.g, p.b), C(p.p.x, p.p.y, p.p.z), rr = right * p.size, uu = up * p.size;
            glm::vec3 v0 = C - rr - uu, v1 = C + rr - uu, v2 = C + rr + uu, v3 = C - rr + uu;
            float qd[6][9] = {{v0.x,v0.y,v0.z,c.r,c.g,c.b,a,0,0},{v1.x,v1.y,v1.z,c.r,c.g,c.b,a,1,0},{v2.x,v2.y,v2.z,c.r,c.g,c.b,a,1,1},
                              {v0.x,v0.y,v0.z,c.r,c.g,c.b,a,0,0},{v2.x,v2.y,v2.z,c.r,c.g,c.b,a,1,1},{v3.x,v3.y,v3.z,c.r,c.g,c.b,a,0,1}};
            for (auto& q : qd) qv.insert(qv.end(), q, q + 9); }
        glBindBuffer(GL_ARRAY_BUFFER, qvbo); glBufferData(GL_ARRAY_BUFFER, qv.size() * sizeof(float), qv.data(), GL_DYNAMIC_DRAW);
        return (GLsizei)(qv.size() / 9);
    };

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 3.0, 0); cam.dist = 5.3f; cam.yaw = 0.32f; cam.pitch = 0.05f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    double t = 0;
    auto emitParts = [&]() {
        for (auto& L : leaves) {
            if (!L.ignited) continue;
            for (int c = 0; c < L.GW * L.GH; c++) {
                if (!L.isBurning(c)) continue;
                int i = c % L.GW, j = c / L.GW; Vector3 wp = L.node(i, j);
                if (rnd() < 0.020f) smoke.push_back({wp, Vector3((rnd() - 0.5f) * 0.18, 0.20 + rnd() * 0.22, (rnd() - 0.5f) * 0.18), 0, 0.7f + rnd() * 0.4f, 0.02f, 0.14f, 0.13f, 0.12f, 0.10f});
                if (rnd() < 0.03f) ember.push_back({wp, Vector3((rnd() - 0.5f) * 0.13, 0.05 + rnd() * 0.13, (rnd() - 0.5f) * 0.13), 0, 0.30f + rnd() * 0.28f, 0.030f, 1.4f, 0.42f, 0.05f, 0.11f});
            }
        }
    };
    auto updateParts = [&](std::vector<Part>& ps, real dt, bool grow) {
        for (auto& p : ps) { p.life += dt; p.v.y += (grow ? 0.15 : 0.2) * dt; p.v *= 0.99; p.p += p.v * dt; if (grow) p.size += 0.10f * dt; }
        ps.erase(std::remove_if(ps.begin(), ps.end(), [](const Part& p) { return p.life >= p.max; }), ps.end());
    };
    auto stepPhysics = [&](real dt) {
        Vector3 wind(std::sin(t * 0.6) * 0.35 + 0.1, 0, std::cos(t * 0.45) * 0.35);
        for (auto& L : leaves) { L.step(dt, wind); if (L.consumed() || L.pos.y < -0.3) spawn(L, false); }
        emitParts(); updateParts(smoke, dt, true); updateParts(ember, dt, false);
        t += dt;
    };

    auto renderFrame = [&]() {
        buildLeaves();
        r.setLightForScene(glm::vec3(0, 2.0, 0), 6.0f);
        r.beginShadow();
        gfx::setM4(r.pDepth, "uModel", glm::mat4(1)); glBindVertexArray(lvao); glDrawArrays(GL_TRIANGLES, 0, leafVerts);
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.19, 0.15, 0.10), 0.95f, 0.0f);
        glDisable(GL_CULL_FACE);
        glUseProgram(pLeaf); gfx::setM4(pLeaf, "uView", cam.view()); gfx::setM4(pLeaf, "uProj", proj); gfx::setV3(pLeaf, "uLightDir", r.lightDir);
        glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D_ARRAY, leafTex); gfx::setI(pLeaf, "uLeaf", 1); glActiveTexture(GL_TEXTURE0);
        glBindVertexArray(lvao); glDrawArrays(GL_TRIANGLES, 0, leafVerts);
        glm::vec3 eye = cam.eye(), fwd = glm::normalize(cam.target - eye);
        glm::vec3 rt = glm::normalize(glm::cross(fwd, glm::vec3(0, 1, 0))), upv = glm::cross(rt, fwd);
        glUseProgram(r.pPart); gfx::setM4(r.pPart, "uView", cam.view()); gfx::setM4(r.pPart, "uProj", proj);
        glEnable(GL_BLEND); glDepthMask(GL_FALSE);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        { GLsizei n = buildQuads(smoke, rt, upv); glBindVertexArray(qvao); glDrawArrays(GL_TRIANGLES, 0, n); }
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        { GLsizei n = buildQuads(ember, rt, upv); glBindVertexArray(qvao); glDrawArrays(GL_TRIANGLES, 0, n); }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 0.35f + 0.1f * std::sin(f * 6.2831853 / frames); renderFrame();
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
