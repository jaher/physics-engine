// A room full of objects and an explosion. Several solid objects — a block, some
// pillars, a slab — sit in an enclosed room, each PRE-FRACTURED with Voronoi cells
// (phys::voronoi) into irregular convex chunks that stay welded until the blast. A
// charge at the centre detonates and shatters every object at once, hurling the
// convex fragments across the room where they ricochet off the six walls (convex
// hull vs half-space, so nothing tunnels out) and settle into rubble. The room is a
// cutaway (near walls hidden); the blast is an emissive fireball with a dust plume.
//   ./roomblast3d --shot out.png [frames]   |   ./roomblast3d --video f [n]
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
struct Frag { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionConvex> hull; gfx::Mesh mesh; float rad; glm::vec3 colour; };

// convex fragment vs the room's 6 wall half-spaces, plus convex-convex among fragments.
struct Collider : ContactGenerator {
    std::vector<Frag>* fr; const CollisionPlane* walls; int nWalls; bool* shattered; real fri, res;
    unsigned addContact(Contact* c, unsigned limit) const override {
        if (!*shattered) return 0;
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fri; d.restitution = res;
        for (auto& f : *fr) f.hull->calculateInternals();
        for (auto& f : *fr) for (int w = 0; w < nWalls; w++) { if (!d.hasMoreContacts()) break; ConvexCollision::convexAndHalfSpace(*f.hull, walls[w], &d); }
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

struct Puff { glm::vec3 p, v; float life, life0, size, grow, buoy; bool fire; };
static float smooth01(float a, float b, float x) { float t = glm::clamp((x - a) / (b - a), 0.0f, 1.0f); return t * t * (3 - 2 * t); }
static glm::vec3 fireColor(float t) {
    glm::vec3 c = glm::mix(glm::vec3(0.30, 0.03, 0.0), glm::vec3(1.05, 0.16, 0.02), smooth01(0.0, 0.35, t));
    c = glm::mix(c, glm::vec3(1.6, 0.62, 0.10), smooth01(0.30, 0.60, t));
    c = glm::mix(c, glm::vec3(2.0, 1.45, 0.45), smooth01(0.55, 0.82, t));
    c = glm::mix(c, glm::vec3(2.6, 2.35, 2.0), smooth01(0.80, 1.0, t));
    return c;
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 360;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    RigidBody::sleepEpsilon = (real)0.7;

    // --- room half-widths and height as 6 inward half-spaces ---
    const real Rx = 6.8, Rz = 6.8, Hy = 5.4;
    CollisionPlane walls[6];
    walls[0].direction = Vector3(0, 1, 0);  walls[0].offset = 0;
    walls[1].direction = Vector3(0, -1, 0); walls[1].offset = -Hy;
    walls[2].direction = Vector3(-1, 0, 0); walls[2].offset = -Rx;
    walls[3].direction = Vector3(1, 0, 0);  walls[3].offset = -Rx;
    walls[4].direction = Vector3(0, 0, -1); walls[4].offset = -Rz;
    walls[5].direction = Vector3(0, 0, 1);  walls[5].offset = -Rz;

    // --- objects, each Voronoi-fractured into welded convex fragments ---
    std::vector<Frag> frags; std::vector<ConvexCell> cells; bool shattered = false;
    auto addObject = [&](Vector3 center, Vector3 half, int nseed, real density, glm::vec3 base, unsigned seed) {
        unsigned s = seed; auto rr = [&]() { s = s * 1103515245u + 12345u; return (float)(((s >> 16) & 0x7fff) / 32767.0); };
        std::vector<Vector3> seeds;
        for (int i = 0; i < nseed; i++) seeds.push_back(Vector3((rr() * 2 - 1) * half.x * 0.92, (rr() * 2 - 1) * half.y * 0.92, (rr() * 2 - 1) * half.z * 0.92));
        for (auto& cell : voronoiFracture(half, seeds, 22)) {
            Frag f; f.rad = 0; for (auto& v : cell.verts) f.rad = std::max(f.rad, (float)v.magnitude());
            f.body = std::make_unique<RigidBody>(); f.body->setMass(cell.volume * density);
            Matrix3 I; for (int k = 0; k < 9; k++) I.data[k] = cell.inertiaUnit.data[k] * density; f.body->setInertiaTensor(I);
            f.body->setPosition(center + cell.site); f.body->setAcceleration(0, -9.81, 0);
            f.body->setDamping(0.5, 0.5); f.body->calculateDerivedData(); f.body->setAwake(false);
            f.hull = std::make_unique<CollisionConvex>(); f.hull->body = f.body.get(); f.hull->vertices = cell.verts;
            float t = (rr() - 0.5f) * 0.09f; f.colour = base + glm::vec3(t, t * 0.9f, t * 0.8f);
            frags.push_back(std::move(f)); cells.push_back(std::move(cell));
        }
    };
    glm::vec3 grey(0.50, 0.49, 0.47), sand(0.60, 0.49, 0.31), brick(0.54, 0.28, 0.20), slate(0.40, 0.41, 0.45);
    addObject(Vector3(0, 0.85, 0),      Vector3(0.85, 0.85, 0.85), 36, 2.4, grey, 7u);     // central block
    addObject(Vector3(2.9, 1.5, 1.7),   Vector3(0.40, 1.50, 0.40), 24, 2.2, sand, 13u);    // pillar
    addObject(Vector3(-2.7, 1.5, -2.0), Vector3(0.40, 1.50, 0.40), 24, 2.2, sand, 21u);    // pillar
    addObject(Vector3(-3.3, 0.65, 2.5), Vector3(0.65, 0.65, 0.65), 26, 2.3, brick, 29u);   // brick block
    addObject(Vector3(3.3, 0.70, -2.7), Vector3(0.70, 0.70, 0.70), 26, 2.4, grey, 33u);    // block
    addObject(Vector3(0, 0.40, 3.6),    Vector3(1.50, 0.40, 0.80), 34, 2.3, slate, 41u);   // slab
    addObject(Vector3(-3.6, 1.3, -0.4), Vector3(0.40, 1.30, 0.40), 22, 2.2, sand, 47u);    // pillar

    World world(20000);
    for (auto& f : frags) world.getRigidBodies().push_back(f.body.get());
    std::vector<char> locked(frags.size(), 0); std::vector<int> pcalm(frags.size(), 0);
    Collider col; col.fr = &frags; col.walls = walls; col.nWalls = 6; col.shattered = &shattered; col.fri = 0.62; col.res = 0.12;
    world.getContactGenerators().push_back(&col);

    // --- blast fire/smoke ---
    std::vector<Puff> puffs;
    unsigned pc = 91u; auto cr = [&]() { pc = pc * 1103515245u + 12345u; return (float)(((pc >> 16) & 0x7fff) / 32767.0); };
    Vector3 charge(0, 0.9, 0); glm::vec3 chargeG(0, 0.9f, 0);
    bool detonated = false; double simT = 0; const double BOOM = 0.5;

    auto detonate = [&]() {
        shattered = true; unsigned st = 5u; auto rr = [&]() { st = st * 1103515245u + 12345u; return (real)(((st >> 16) & 0x7fff) / 32767.0); };
        for (auto& f : frags) {                                       // fierce radial blast that shatters every object
            Vector3 dr = f.body->getPosition() - charge; real dist = dr.magnitude();
            Vector3 dir = dist > 1e-6 ? dr * (1.0 / dist) : Vector3(0, 1, 0);
            real speed = 15.0 / (dist + 0.5);
            Vector3 v = dir * speed + Vector3(0, 0.5 * speed + 2.0, 0) + Vector3(rr() - 0.5, rr() - 0.5, rr() - 0.5) * 2.5;
            f.body->setAwake(true); f.body->setVelocity(v);
            f.body->setRotation(Vector3(rr() - 0.5, rr() - 0.5, rr() - 0.5) * 18.0);
        }
        auto rdir = [&]() { glm::vec3 d(cr() - 0.5f, cr() - 0.5f, cr() - 0.5f); float m = glm::length(d); return m > 1e-4f ? d / m : glm::vec3(0, 1, 0); };
        for (int i = 0; i < 230; i++) {                               // fireball
            Puff q; q.fire = true; q.p = chargeG + rdir() * (0.2f + cr() * 1.5f);
            q.v = rdir() * (4.0f + cr() * 11.0f) + glm::vec3(0, 2.5f + cr() * 3.5f, 0);
            q.life = q.life0 = 0.32f + cr() * 0.7f; q.size = 0.22f + cr() * 0.55f; q.grow = 1.8f + cr() * 1.6f; q.buoy = 4.6f + cr() * 3.0f;
            puffs.push_back(q);
        }
        for (int i = 0; i < 210; i++) {                               // dust plume
            Puff q; q.fire = false; q.p = chargeG + rdir() * (0.15f + cr() * 1.3f);
            q.v = rdir() * (1.4f + cr() * 3.5f) + glm::vec3(0, 2.4f + cr() * 2.0f, 0);
            q.life = q.life0 = 1.6f + cr() * 2.4f; q.size = 0.24f + cr() * 0.45f; q.grow = 0.9f + cr() * 0.7f; q.buoy = 2.0f + cr() * 1.4f;
            puffs.push_back(q);
        }
    };
    auto stepPhysics = [&]() {
        for (int k = 0; k < 2; k++) { world.startFrame(); world.runPhysics(1.0 / 240); }
        simT += 1.0 / 120;
        if (!detonated && simT >= BOOM) { detonate(); detonated = true; }
        else if (detonated) {
            double td = simT - BOOM;                                  // bleed slow fragments, then lock them still
            for (size_t i = 0; i < frags.size(); i++) {
                RigidBody* b = frags[i].body.get();
                if (locked[i]) { b->setVelocity(Vector3()); b->setRotation(Vector3()); continue; }
                real v = b->getVelocity().magnitude(), w = b->getRotation().magnitude();
                if (v < 0.6) { b->setVelocity(b->getVelocity() * 0.82); b->setRotation(b->getRotation() * 0.70); }
                if (v < 0.35 && w < 0.9) pcalm[i]++; else pcalm[i] = 0;
                if (pcalm[i] > 14 || td > 4.4) {
                    b->setVelocity(Vector3()); b->setRotation(Vector3()); b->setAcceleration(0, 0, 0);
                    b->setInverseMass(0); Matrix3 bi; bi.setInertiaTensorCoeffs(1e12, 1e12, 1e12); b->setInertiaTensor(bi);
                    b->calculateDerivedData(); b->setAwake(false); locked[i] = 1;
                }
            }
        }
        float dt = 1.0f / 120;
        for (auto& q : puffs) { q.v.y += q.buoy * dt; q.v *= (q.fire ? 0.90f : 0.95f); q.p += q.v * dt; q.size += q.grow * dt; q.life -= dt; }
        puffs.erase(std::remove_if(puffs.begin(), puffs.end(), [](const Puff& q) { return q.life <= 0; }), puffs.end());
    };

    gfx::App app(W, H, "room explosion", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), plane = gfx::makePlane(20, 1);
    for (size_t i = 0; i < frags.size(); i++) {                       // build each fragment's flat-shaded mesh (needs GL)
        auto& cell = cells[i]; std::vector<float> vv; std::vector<unsigned> ii;
        for (size_t t = 0; t < cell.tris.size(); t++) { Vector3 n = cell.triN[t];
            for (int k = 0; k < 3; k++) { Vector3 p = cell.verts[cell.tris[t][k]];
                vv.insert(vv.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z}); ii.push_back((unsigned)ii.size()); } }
        frags[i].mesh.upload(vv, ii);
    }
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.4, 0); cam.dist = 15.5f; cam.yaw = 0.6f; cam.pitch = 0.32f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    // billboard fire/smoke program (drawn into the renderer's HDR MSAA buffer)
    GLuint pBB = gfx::program(R"(#version 330 core
layout(location=0) in vec3 aC; layout(location=1) in vec2 aUV; layout(location=2) in float aS; layout(location=3) in vec4 aCol;
uniform mat4 uV,uP; uniform vec3 uRight,uUp; out vec2 vUV; out vec4 vCol; out vec2 vSeed;
void main(){ vUV=aUV; vCol=aCol; vSeed=aC.xz; vec3 wp=aC+(aUV.x*uRight+aUV.y*uUp)*aS; gl_Position=uP*uV*vec4(wp,1.0); })",
        R"(#version 330 core
in vec2 vUV; in vec4 vCol; in vec2 vSeed; out vec4 frag; uniform int uAdd;
float hash(vec2 p){ return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453); }
float vnoise(vec2 p){ vec2 i=floor(p),f=fract(p); f=f*f*(3.0-2.0*f);
  float a=hash(i),b=hash(i+vec2(1,0)),c=hash(i+vec2(0,1)),d=hash(i+vec2(1,1));
  return mix(mix(a,b,f.x),mix(c,d,f.x),f.y); }
void main(){ float rr=length(vUV); if(rr>1.0) discard;
  float n = vnoise(vUV*2.6+vSeed*4.1)*0.65 + vnoise(vUV*6.0+vSeed*2.3)*0.35;
  float soft = smoothstep(1.0,0.12,rr) * (0.12+1.15*n); soft = clamp(soft,0.0,1.0);
  if(soft<0.06) discard;
  if(uAdd==1) frag=vec4(vCol.rgb*soft*vCol.a, 1.0);
  else        frag=vec4(vCol.rgb, soft*vCol.a); })");
    GLuint bbVAO, bbVBO; glGenVertexArrays(1, &bbVAO); glGenBuffers(1, &bbVBO);
    glBindVertexArray(bbVAO); glBindBuffer(GL_ARRAY_BUFFER, bbVBO);
    glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2); glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)(5 * sizeof(float)));
    glEnableVertexAttribArray(3); glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)(6 * sizeof(float)));
    glBindVertexArray(0);

    auto emitBatch = [&](std::vector<float>& v, const Puff& q) {
        float t = q.life / q.life0; glm::vec4 col;
        if (q.fire) { float temp = glm::clamp(t * 1.15f, 0.0f, 1.0f); glm::vec3 c = fireColor(temp) * (0.35f + 0.85f * temp); col = glm::vec4(c, glm::clamp(t * 1.3f, 0.0f, 1.0f)); }
        else { float g = 0.16f + 0.10f * t; float shade = smooth01(0.0f, 0.4f, t); col = glm::vec4(g, g * 0.96f, g * 0.90f, 0.34f * shade); }
        const float uv[6][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, -1}, {1, 1}, {-1, 1}};
        for (auto& c : uv) v.insert(v.end(), {q.p.x, q.p.y, q.p.z, c[0], c[1], q.size, col.r, col.g, col.b, col.a});
    };
    auto drawPuffs = [&](glm::vec3 eye, glm::vec3 R, glm::vec3 U) {
        std::vector<float> vs, vf; glm::vec3 fwd = glm::normalize(cam.target - eye);
        std::vector<const Puff*> sm; for (auto& q : puffs) if (!q.fire) sm.push_back(&q);
        std::sort(sm.begin(), sm.end(), [&](const Puff* a, const Puff* b) { return glm::dot(a->p - eye, fwd) > glm::dot(b->p - eye, fwd); });
        for (auto* q : sm) emitBatch(vs, *q);
        for (auto& q : puffs) if (q.fire) emitBatch(vf, q);
        glUseProgram(pBB);
        glUniformMatrix4fv(glGetUniformLocation(pBB, "uV"), 1, GL_FALSE, &r.curView[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(pBB, "uP"), 1, GL_FALSE, &r.curProj[0][0]);
        glUniform3fv(glGetUniformLocation(pBB, "uRight"), 1, &R[0]); glUniform3fv(glGetUniformLocation(pBB, "uUp"), 1, &U[0]);
        glBindVertexArray(bbVAO); glEnable(GL_BLEND); glDepthMask(GL_FALSE);
        if (!vs.empty()) { glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); glUniform1i(glGetUniformLocation(pBB, "uAdd"), 0);
            glBindBuffer(GL_ARRAY_BUFFER, bbVBO); glBufferData(GL_ARRAY_BUFFER, vs.size() * sizeof(float), vs.data(), GL_DYNAMIC_DRAW); glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(vs.size() / 10)); }
        if (!vf.empty()) { glBlendFunc(GL_ONE, GL_ONE); glUniform1i(glGetUniformLocation(pBB, "uAdd"), 1);
            glBindBuffer(GL_ARRAY_BUFFER, bbVBO); glBufferData(GL_ARRAY_BUFFER, vf.size() * sizeof(float), vf.data(), GL_DYNAMIC_DRAW); glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(vf.size() / 10)); }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
    };

    struct WallDraw { glm::vec3 c, s, nIn, col; };
    std::vector<WallDraw> wd = {
        {{Rx + 0.1f, Hy * 0.5f, 0}, {0.1f, Hy * 0.5f, Rz}, {-1, 0, 0}, {0.52, 0.52, 0.55}},
        {{-Rx - 0.1f, Hy * 0.5f, 0}, {0.1f, Hy * 0.5f, Rz}, {1, 0, 0}, {0.52, 0.52, 0.55}},
        {{0, Hy * 0.5f, Rz + 0.1f}, {Rx, Hy * 0.5f, 0.1f}, {0, 0, -1}, {0.55, 0.54, 0.52}},
        {{0, Hy * 0.5f, -Rz - 0.1f}, {Rx, Hy * 0.5f, 0.1f}, {0, 0, 1}, {0.55, 0.54, 0.52}},
        {{0, Hy + 0.1f, 0}, {Rx, 0.1f, Rz}, {0, -1, 0}, {0.42, 0.42, 0.46}},
    };

    auto renderFrame = [&]() {
        glm::vec3 eye = cam.eye(); glm::vec3 fwd = glm::normalize(cam.target - eye);
        glm::vec3 R = glm::normalize(glm::cross(fwd, glm::vec3(0, 1, 0))), U = glm::cross(R, fwd);
        r.setLightForScene(glm::vec3(0, 1.6, 0), 11.0f);
        r.beginShadow();
        for (auto& f : frags) r.shadowDraw(f.mesh, glmFromPhys(f.body->getTransform()));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        glClearColor(0.15f, 0.16f, 0.19f, 1.0f); glClear(GL_COLOR_BUFFER_BIT);
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.30, 0.29, 0.28), 0.96f, 0.0f);
        for (auto& w : wd) if (glm::dot(eye - w.c, w.nIn) > 0)
            r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), w.c), w.s), w.col, 0.92f, 0.0f);
        for (auto& f : frags) r.drawPBR(r.pSolid, f.mesh, glmFromPhys(f.body->getTransform()), f.colour, 0.8f, 0.0f);
        drawPuffs(eye, R, U);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw += 0.5f / frames; renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 2; s++) stepPhysics(); }
        std::printf("wrote %d frames (%zu fragments)\n", frames, frags.size()); return 0;
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
