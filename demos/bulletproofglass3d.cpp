// A bullet hitting bulletproof (laminated) glass. Unlike a plain pane, ballistic
// glass STOPS the round and holds together — the polymer interlayer bonds the
// fragments so nothing flies off. The realistic signature is a localized milky-white
// crater: a pulverized star at the impact ringed by radial + concentric cracks that
// fade into clear glass, the pane staying whole in its frame while the spent, flattened
// bullet drops down the face. We reuse the impact fracture pattern (phys::glassfrac)
// for the crack web, but the shards never detach — only their look changes, a white
// craze spreading outward from the hit — and the bullet is arrested at the surface.
//   ./bulletproofglass3d --shot out.png [frames]   |   ./bulletproofglass3d --video f [n]
#include "phys/phys.h"
#include "common/gfx.h"
#include <vector>
#include <cstring>
#include <cmath>
#include <algorithm>
using namespace phys;

struct Shard { gfx::Mesh mesh; glm::vec3 pos; float dist; };
struct Dust { glm::vec3 p, v; float life, life0, size; };

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 240;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    // glass pane geometry & impact point (upper-centre so the spent bullet drops in view)
    const real Cy = 1.45, PW = 1.5, PH = 1.0, PT = 0.05;             // thicker → ballistic laminate
    const real iu = 0.10, iv = 0.32;
    glm::vec3 impact((float)iu, (float)(Cy + iv), 0);
    const float glassFront = (float)(PT * 0.5f);

    // pre-fracture for the crack web; the shards are STATIC (never detach)
    auto cells = glassFracture(PW, PH, PT, iu, iv, 30, 10, 0.035, 7u);
    std::vector<Shard> shards;
    for (auto& cell : cells) {
        Shard s; s.pos = glm::vec3((float)cell.site.x, (float)(Cy + cell.site.y), (float)cell.site.z);
        s.dist = (float)std::hypot((double)(cell.site.x - iu), (double)(cell.site.y - iv));
        shards.push_back(std::move(s));
    }

    // bullet (kinematic): flies in, is arrested at the glass, then drops down the face
    glm::vec3 bpos(impact.x, impact.y, 3.0f); glm::vec3 bvel(0, 0, -15.0f); const float bR = 0.055f;
    bool stopped = false; double simT = 0, impactT = -1; float bflat = 1.0f;
    std::vector<Dust> dust;
    unsigned ds = 5u; auto rr = [&]() { ds = ds * 1103515245u + 12345u; return (float)(((ds >> 16) & 0x7fff) / 32767.0f); };

    auto stepPhysics = [&]() {
        double dt = 1.0 / 120; simT += dt;
        if (!stopped) {
            bpos += bvel * (float)dt;
            if (bpos.z <= glassFront) {                              // arrested by the ballistic glass
                stopped = true; impactT = simT; bpos.z = glassFront + 0.02f; bvel = glm::vec3(0.15f, -0.4f, 0);
                bflat = 0.4f;                                        // mushroomed round
                for (int i = 0; i < 60; i++) {                       // pulverized-glass spray off the face
                    Dust d; d.p = impact + glm::vec3(0, 0, glassFront + 0.02f);
                    glm::vec3 dir(rr() - 0.5f, rr() - 0.5f, -(0.3f + rr()));      // mostly toward the shooter
                    d.v = glm::normalize(dir) * (0.6f + rr() * 2.2f);
                    d.life = d.life0 = 0.3f + rr() * 0.5f; d.size = 0.012f + rr() * 0.02f;
                    dust.push_back(d);
                }
            }
        } else {                                                    // spent bullet slides/drops down the glass
            bvel.y -= 6.5f * (float)dt; bvel.x *= 0.98f; bpos += bvel * (float)dt;
            if (bpos.y < bR) { bpos.y = bR; bvel = glm::vec3(0); }
        }
        for (auto& d : dust) { d.v.y -= 5.0f * (float)dt; d.p += d.v * (float)dt; d.life -= (float)dt; }
        dust.erase(std::remove_if(dust.begin(), dust.end(), [](const Dust& d) { return d.life <= 0; }), dust.end());
    };
    // per-shard crazing: 0 (clear) → 1 (milky white), spreading outward from the hit
    const float crazedR = 0.62f;
    auto damageOf = [&](const Shard& s) -> float {
        if (impactT < 0) return 0.0f;
        float target = std::pow(std::max(0.0f, 1.0f - s.dist / crazedR), 0.55f);
        float te = (float)(simT - impactT) - s.dist * 0.10f;         // shockwave reaches farther shards later
        float grow = glm::clamp(te / 0.09f, 0.0f, 1.0f); grow = grow * grow * (3 - 2 * grow);
        return target * grow;
    };

    gfx::App app(W, H, "bulletproof glass", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(30, 1);
    for (size_t i = 0; i < shards.size(); i++) {
        auto& cell = cells[i]; std::vector<float> vv; std::vector<unsigned> ii;
        for (size_t t = 0; t < cell.tris.size(); t++) { Vector3 n = cell.triN[t];
            for (int k = 0; k < 3; k++) { Vector3 p = cell.verts[cell.tris[t][k]];
                vv.insert(vv.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z}); ii.push_back((unsigned)ii.size()); } }
        shards[i].mesh.upload(vv, ii);
    }
    gfx::Mesh pane; { std::vector<float> v = {-PW, -PH, 0, 0, 0, -1, PW, -PH, 0, 0, 0, -1, PW, PH, 0, 0, 0, -1, -PW, PH, 0, 0, 0, -1};
        std::vector<unsigned> idx = {0, 1, 2, 0, 2, 3}; pane.upload(v, idx); }

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.4, 0); cam.dist = 5.2f; cam.yaw = 2.55f; cam.pitch = 0.11f;
    glm::mat4 proj = glm::perspective(glm::radians(43.0f), (float)W / H, 0.05f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    // glass shader with a per-draw damage (clear → milky-white craze, bright cracks)
    GLuint pGlass = gfx::program(R"(#version 330 core
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aN;
uniform mat4 uM,uV,uP; out vec3 vW; out vec3 vN;
void main(){ vec4 w=uM*vec4(aPos,1.0); vW=w.xyz; vN=mat3(uM)*aN; gl_Position=uP*uV*w; })",
        R"(#version 330 core
in vec3 vW; in vec3 vN; out vec4 frag; uniform vec3 uCam,uLight; uniform float uDmg;
void main(){ vec3 N=normalize(vN); vec3 V=normalize(uCam-vW); if(dot(N,V)<0.0) N=-N;
  vec3 L=normalize(uLight);
  float fres=pow(1.0-abs(dot(N,V)),3.0);
  float spec=pow(max(dot(reflect(-L,N),V),0.0),80.0);
  vec3 tint=mix(vec3(0.55,0.73,0.82), vec3(0.92,0.94,0.97), uDmg);   // clear → milky white
  float diff=max(dot(N,L),0.0);
  vec3 crack=vec3(0.9,0.95,1.0)*fres*(0.55 + 2.2*uDmg);              // seams glow white where crazed
  vec3 col=tint*(0.32+0.4*diff) + crack + vec3(spec)*1.1;
  float a=clamp(mix(0.10,0.94,uDmg) + fres*(0.35+0.4*uDmg) + spec, 0.0, 0.97);
  frag=vec4(col,a); })");
    auto setM = [](GLuint p, const char* n, const glm::mat4& m) { glUniformMatrix4fv(glGetUniformLocation(p, n), 1, GL_FALSE, &m[0][0]); };
    auto setV = [](GLuint p, const char* n, const glm::vec3& v) { glUniform3fv(glGetUniformLocation(p, n), 1, &v[0]); };
    auto setF = [](GLuint p, const char* n, float f) { glUniform1f(glGetUniformLocation(p, n), f); };

    auto renderFrame = [&]() {
        glm::vec3 eye = cam.eye();
        r.setLightForScene(glm::vec3(0, 1.3, 0.5), 6.0f);
        r.beginShadow();
        r.shadowDraw(sphere, glm::scale(glm::translate(glm::mat4(1), bpos), glm::vec3(bR)));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        glClearColor(0.20f, 0.23f, 0.28f, 1.0f); glClear(GL_COLOR_BUFFER_BIT);
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.34, 0.33, 0.32), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, 1.6, 3.4)), glm::vec3(4.5, 2.2, 0.1)), glm::vec3(0.40, 0.36, 0.30), 0.9f, 0.0f);
        glm::vec3 fc(0.12, 0.10, 0.09); float fb = 0.10f;            // window frame
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, Cy + PH + fb, 0)), glm::vec3(PW + fb, fb, 0.10)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, Cy - PH - fb, 0)), glm::vec3(PW + fb, fb, 0.10)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(-PW - fb, Cy, 0)), glm::vec3(fb, PH, 0.10)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(PW + fb, Cy, 0)), glm::vec3(fb, PH, 0.10)), fc, 0.7f, 0.1f);
        glm::mat4 bm = glm::scale(glm::translate(glm::mat4(1), bpos), glm::vec3(bR, bR, bR * bflat));   // flattened after impact
        r.drawPBR(r.pSolid, sphere, bm, glm::vec3(0.05, 0.05, 0.06), 0.3f, 0.8f);
        for (auto& d : dust) { float a = d.life / d.life0; r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), d.p), glm::vec3(d.size * (0.4f + a))), glm::vec3(0.85, 0.9, 0.95), 0.6f, 0.0f); }

        // --- glass pass (transparent; static shards, crazed near the hit) ---
        glUseProgram(pGlass);
        setM(pGlass, "uV", r.curView); setM(pGlass, "uP", r.curProj);
        setV(pGlass, "uCam", eye); setV(pGlass, "uLight", glm::normalize(glm::vec3(-0.4, 0.7, 0.5)));
        glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); glDepthMask(GL_FALSE);
        if (impactT < 0) { setF(pGlass, "uDmg", 0.0f); setM(pGlass, "uM", glm::translate(glm::mat4(1), glm::vec3(0, Cy, 0))); pane.draw(); }
        else {
            std::vector<int> order(shards.size()); for (size_t i = 0; i < shards.size(); i++) order[i] = (int)i;
            std::sort(order.begin(), order.end(), [&](int a, int b) { return glm::length(shards[a].pos - eye) > glm::length(shards[b].pos - eye); });
            for (int i : order) { setF(pGlass, "uDmg", damageOf(shards[i])); setM(pGlass, "uM", glm::translate(glm::mat4(1), shards[i].pos)); shards[i].mesh.draw(); }
        }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw += 0.3f / frames; renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 2; s++) stepPhysics(); }
        std::printf("wrote %d frames (%zu shards, held)\n", frames, shards.size()); return 0;
    }
    if (headless) { for (int i = 0; i < frames * 2; i++) stepPhysics(); renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        stepPhysics(); renderFrame(); r.present(); app.poll();
    }
    return 0;
}
