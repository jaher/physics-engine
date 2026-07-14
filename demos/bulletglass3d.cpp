// A bullet through a pane of glass. The pane is pre-fractured with an impact
// ("bullet-hole") pattern (phys::glassfrac) centred exactly where the bullet will
// strike — a punched hole, radial spokes, concentric rings, fine shards at the entry
// and coarse plates at the rim — but the shards stay welded so the window looks
// whole. When the bullet reaches the glass they detach: the shards around the hole
// are ejected forward in a spall cone while the rest cracks free and falls, all as
// thin convex fragments that tumble and land. The glass is drawn as a real
// transparent, Fresnel-edged material, depth-sorted, over the room behind it.
//   ./bulletglass3d --shot out.png [frames]   |   ./bulletglass3d --video f [n]
#include "phys/phys.h"
#include "common/gfx.h"
#include <vector>
#include <memory>
#include <cstring>
#include <cmath>
#include <algorithm>
using namespace phys;

static glm::mat4 glmFromPhys(const Matrix4& t) {
    return glm::mat4(t.data[0], t.data[4], t.data[8], 0, t.data[1], t.data[5], t.data[9], 0,
                     t.data[2], t.data[6], t.data[10], 0, t.data[3], t.data[7], t.data[11], 1);
}
struct Frag { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionConvex> hull; gfx::Mesh mesh; float rad; glm::vec3 c0; };

struct Collider : ContactGenerator {
    std::vector<Frag>* fr; CollisionPlane floor; bool* shattered; real fri, res;
    unsigned addContact(Contact* c, unsigned limit) const override {
        if (!*shattered) return 0;
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fri; d.restitution = res;
        for (auto& f : *fr) f.hull->calculateInternals();
        for (auto& f : *fr) { if (!d.hasMoreContacts()) break; ConvexCollision::convexAndHalfSpace(*f.hull, floor, &d); }
        for (size_t i = 0; i < fr->size(); i++) { if (!d.hasMoreContacts()) break;
            Vector3 pi = (*fr)[i].body->getPosition(); float ri = (*fr)[i].rad;
            for (size_t j = i + 1; j < fr->size(); j++) { if (!d.hasMoreContacts()) break;
                Vector3 pj = (*fr)[j].body->getPosition(); float rj = (*fr)[j].rad;
                if ((pi - pj).squareMagnitude() > (ri + rj) * (ri + rj)) continue;
                ConvexCollision::convexAndConvex(*(*fr)[i].hull, *(*fr)[j].hull, &d);
            }
        }
        return d.contactCount;
    }
};

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    RigidBody::sleepEpsilon = (real)0.4;

    // --- glass pane geometry & impact point ---
    const real Cy = 1.45, PW = 1.5, PH = 1.0, PT = 0.03;             // pane centre-y, half-w, half-h, thickness
    const real iu = 0.12, iv = 0.18;                                 // impact in pane-local (u,v)
    glm::vec3 impact((float)iu, (float)(Cy + iv), 0);

    std::vector<Frag> frags; std::vector<ConvexCell> cells; bool shattered = false;
    auto glassCells = glassFracture(PW, PH, PT, iu, iv, 26, 8, 0.055, 7u);
    const real gdensity = 2.5;
    for (auto& cell : glassCells) {
        Frag f; f.rad = 0; for (auto& v : cell.verts) f.rad = std::max(f.rad, (float)v.magnitude());
        f.body = std::make_unique<RigidBody>(); f.body->setMass(std::max(1e-3, cell.volume * gdensity));
        Matrix3 I; for (int k = 0; k < 9; k++) I.data[k] = cell.inertiaUnit.data[k] * gdensity;
        for (int k : {0, 4, 8}) I.data[k] = std::max(I.data[k], (real)1e-6); f.body->setInertiaTensor(I);
        f.body->setPosition(Vector3(cell.site.x, Cy + cell.site.y, cell.site.z));
        f.body->setAcceleration(0, -9.81, 0); f.body->setDamping(0.6, 0.4); f.body->calculateDerivedData(); f.body->setAwake(false);
        f.hull = std::make_unique<CollisionConvex>(); f.hull->body = f.body.get(); f.hull->vertices = cell.verts;
        frags.push_back(std::move(f)); cells.push_back(std::move(cell));
    }

    World world(16000);
    for (auto& f : frags) world.getRigidBodies().push_back(f.body.get());
    std::vector<char> locked(frags.size(), 0); std::vector<int> pcalm(frags.size(), 0);
    Collider col; col.fr = &frags; col.floor.direction = Vector3(0, 1, 0); col.floor.offset = 0; col.shattered = &shattered; col.fri = 0.5; col.res = 0.1;
    world.getContactGenerators().push_back(&col);

    // --- bullet (kinematic: flies straight along −z, punches through) ---
    glm::vec3 bpos(impact.x, impact.y, 3.0f); glm::vec3 bvel(0, 0, -15.0f); const float bR = 0.05f;
    double simT = 0; bool fired = false;

    auto shatter = [&]() {
        shattered = true; unsigned st = 3u; auto rr = [&]() { st = st * 1103515245u + 12345u; return (real)(((st >> 16) & 0x7fff) / 32767.0); };
        Vector3 fwd(0, 0, -1);                                       // bullet direction
        for (auto& f : frags) {
            Vector3 p = f.body->getPosition();
            real du = p.x - impact.x, dv = p.y - impact.y; real dist = std::sqrt(du * du + dv * dv);
            Vector3 rad = dist > 1e-4 ? Vector3(du, dv, 0) * (1.0 / dist) : Vector3(0, 1, 0);
            real fwdSpeed = 2.6 / (dist + 0.18);                     // spall cone: fast near the hole
            real radSpeed = 1.1 / (dist + 0.25);
            Vector3 v = fwd * fwdSpeed + rad * radSpeed + Vector3(rr() - 0.5, rr() - 0.5, rr() - 0.5) * 0.6;
            f.body->setAwake(true); f.body->setVelocity(v);
            f.body->setRotation(Vector3(rr() - 0.5, rr() - 0.5, rr() - 0.5) * (4.0 + 10.0 / (dist + 0.3)));
        }
    };
    auto stepPhysics = [&]() {
        for (int k = 0; k < 2; k++) { world.startFrame(); world.runPhysics(1.0 / 240); }
        double dt = 1.0 / 120; simT += dt;
        if (!fired) { bpos += bvel * (float)dt; bvel.y -= 1.2f * (float)dt;   // slight drop
            if (bpos.z <= 0.0f) { shatter(); fired = true; } }
        else { bpos += bvel * (float)dt; bvel.y -= 6.0f * (float)dt; }        // bullet continues, arcs down
        if (fired) for (size_t i = 0; i < frags.size(); i++) {               // settle: bleed slow shards, lock at rest
            RigidBody* b = frags[i].body.get();
            if (locked[i]) { b->setVelocity(Vector3()); b->setRotation(Vector3()); continue; }
            real v = b->getVelocity().magnitude(), w = b->getRotation().magnitude();
            if (v < 0.5) { b->setVelocity(b->getVelocity() * 0.82); b->setRotation(b->getRotation() * 0.72); }
            if (v < 0.3 && w < 0.8) pcalm[i]++; else pcalm[i] = 0;
            if (pcalm[i] > 14 || simT > 4.2) { b->setVelocity(Vector3()); b->setRotation(Vector3()); b->setAcceleration(0, 0, 0);
                b->setInverseMass(0); Matrix3 bi; bi.setInertiaTensorCoeffs(1e12, 1e12, 1e12); b->setInertiaTensor(bi);
                b->calculateDerivedData(); b->setAwake(false); locked[i] = 1; }
        }
    };

    gfx::App app(W, H, "bullet through glass", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(30, 1);
    for (size_t i = 0; i < frags.size(); i++) {                       // build each shard's mesh
        auto& cell = cells[i]; std::vector<float> vv; std::vector<unsigned> ii;
        for (size_t t = 0; t < cell.tris.size(); t++) { Vector3 n = cell.triN[t];
            for (int k = 0; k < 3; k++) { Vector3 p = cell.verts[cell.tris[t][k]];
                vv.insert(vv.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z}); ii.push_back((unsigned)ii.size()); } }
        frags[i].mesh.upload(vv, ii);
    }
    // intact pane (a flat quad) shown until the shot lands
    gfx::Mesh pane; { std::vector<float> v = {
        -PW, -PH, 0, 0, 0, -1,  PW, -PH, 0, 0, 0, -1,  PW, PH, 0, 0, 0, -1,  -PW, PH, 0, 0, 0, -1};
        std::vector<unsigned> idx = {0, 1, 2, 0, 2, 3}; pane.upload(v, idx); }

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.35, 0); cam.dist = 5.6f; cam.yaw = 2.5f; cam.pitch = 0.13f;
    glm::mat4 proj = glm::perspective(glm::radians(43.0f), (float)W / H, 0.05f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    // transparent glass shader (drawn into the renderer's HDR buffer after the opaque scene)
    GLuint pGlass = gfx::program(R"(#version 330 core
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aN;
uniform mat4 uM,uV,uP; out vec3 vW; out vec3 vN;
void main(){ vec4 w=uM*vec4(aPos,1.0); vW=w.xyz; vN=mat3(uM)*aN; gl_Position=uP*uV*w; })",
        R"(#version 330 core
in vec3 vW; in vec3 vN; out vec4 frag; uniform vec3 uCam,uLight;
void main(){ vec3 N=normalize(vN); vec3 V=normalize(uCam-vW); if(dot(N,V)<0.0) N=-N;
  vec3 L=normalize(uLight);
  float fres=pow(1.0-abs(dot(N,V)),3.0);
  float spec=pow(max(dot(reflect(-L,N),V),0.0),80.0);
  vec3 tint=vec3(0.55,0.73,0.82);
  vec3 col=tint*(0.35+0.35*max(dot(N,L),0.0)) + vec3(0.8,0.9,1.0)*fres*0.9 + vec3(spec)*1.2;
  float a=clamp(0.12 + 0.82*fres + spec, 0.0, 0.92);
  frag=vec4(col,a); })");
    auto setM = [](GLuint p, const char* n, const glm::mat4& m) { glUniformMatrix4fv(glGetUniformLocation(p, n), 1, GL_FALSE, &m[0][0]); };
    auto setV = [](GLuint p, const char* n, const glm::vec3& v) { glUniform3fv(glGetUniformLocation(p, n), 1, &v[0]); };

    auto renderFrame = [&]() {
        glm::vec3 eye = cam.eye();
        r.setLightForScene(glm::vec3(0, 1.3, 0.5), 6.0f);
        r.beginShadow();
        if (shattered) for (auto& f : frags) r.shadowDraw(f.mesh, glmFromPhys(f.body->getTransform()));
        r.shadowDraw(sphere, glm::scale(glm::translate(glm::mat4(1), bpos), glm::vec3(bR)));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        glClearColor(0.20f, 0.23f, 0.28f, 1.0f); glClear(GL_COLOR_BUFFER_BIT);
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.34, 0.33, 0.32), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, 1.6, 3.4)), glm::vec3(4.5, 2.2, 0.1)), glm::vec3(0.40, 0.36, 0.30), 0.9f, 0.0f);   // back wall
        // window frame around the pane
        glm::vec3 fc(0.12, 0.10, 0.09); float fb = 0.09f;
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, Cy + PH + fb, 0)), glm::vec3(PW + fb, fb, 0.08)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, Cy - PH - fb, 0)), glm::vec3(PW + fb, fb, 0.08)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(-PW - fb, Cy, 0)), glm::vec3(fb, PH, 0.08)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(PW + fb, Cy, 0)), glm::vec3(fb, PH, 0.08)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), bpos), glm::vec3(bR)), glm::vec3(0.05, 0.05, 0.06), 0.3f, 0.8f);   // bullet

        // --- glass pass (transparent, into the still-bound HDR buffer) ---
        glUseProgram(pGlass);
        setM(pGlass, "uV", r.curView); setM(pGlass, "uP", r.curProj);
        setV(pGlass, "uCam", eye); setV(pGlass, "uLight", glm::normalize(glm::vec3(-0.4, 0.7, 0.5)));
        glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); glDepthMask(GL_FALSE);
        if (!shattered) { setM(pGlass, "uM", glm::translate(glm::mat4(1), glm::vec3(0, Cy, 0))); pane.draw(); }
        else {
            std::vector<int> order(frags.size()); for (size_t i = 0; i < frags.size(); i++) order[i] = (int)i;
            std::sort(order.begin(), order.end(), [&](int a, int b) {                 // back-to-front
                return glm::length(glm::vec3(frags[a].body->getPosition().x, frags[a].body->getPosition().y, frags[a].body->getPosition().z) - eye)
                     > glm::length(glm::vec3(frags[b].body->getPosition().x, frags[b].body->getPosition().y, frags[b].body->getPosition().z) - eye); });
            for (int i : order) { setM(pGlass, "uM", glmFromPhys(frags[i].body->getTransform())); frags[i].mesh.draw(); }
        }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw += 0.32f / frames; renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 2; s++) stepPhysics(); }
        std::printf("wrote %d frames (%zu shards)\n", frames, frags.size()); return 0;
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
