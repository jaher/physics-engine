// Explosion of a concrete block. A solid concrete cube is pre-split into welded
// box fragments (phys::fracture); a charge at its core detonates, throwing every
// fragment radially outward with a blast-overpressure falloff and an upward plume
// (Destructible::detonate). From then on the chunks are ordinary rigid bodies —
// they tumble, collide, and settle into a rubble pile via the engine's contact
// solver. The blast is a billboard fire + smoke system: a white-hot core cools
// through orange to red and lofts a dark smoke plume, drawn additively into the
// renderer's HDR buffer so the tonemapper blooms the core.
//   ./explosion3d --shot out.png [frames]   |   ./explosion3d --video f [n]
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

// convex-convex contacts (GJK/EPA) against the ground box and each other fragment,
// with a broad-phase distance cull so the irregular hulls stay affordable.
struct Collider : ContactGenerator {
    std::vector<Frag>* frags; CollisionConvex* ground; Destructible* dest; real fr, re; const bool* frozen = nullptr;
    unsigned addContact(Contact* c, unsigned limit) const override {
        if (!dest->shattered || (frozen && *frozen)) return 0;
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
        for (auto& f : *frags) f.hull->calculateInternals();
        ground->calculateInternals();
        for (auto& f : *frags) { if (!d.hasMoreContacts()) break; ConvexCollision::convexAndConvex(*f.hull, *ground, &d); }
        for (size_t i = 0; i < frags->size(); i++) { if (!d.hasMoreContacts()) break;
            Vector3 pi = (*frags)[i].body->getPosition(); float ri = (*frags)[i].rad;
            for (size_t j = i + 1; j < frags->size(); j++) { if (!d.hasMoreContacts()) break;
                Vector3 pj = (*frags)[j].body->getPosition(); float rj = (*frags)[j].rad;
                if ((pi - pj).squareMagnitude() > (ri + rj) * (ri + rj)) continue;   // broad-phase cull
                ConvexCollision::convexAndConvex(*(*frags)[i].hull, *(*frags)[j].hull, &d);
            }
        }
        return d.contactCount;
    }
};

// fire/smoke puff. `temp`∈[0,1] cools a fire puff; smoke carries a grey alpha.
struct Puff { glm::vec3 p, v; float life, life0, size, grow, buoy; bool fire; };

static float smooth01(float a, float b, float x) { float t = glm::clamp((x - a) / (b - a), 0.0f, 1.0f); return t * t * (3 - 2 * t); }
static glm::vec3 fireColor(float t) {                        // blackbody-ish: ember → red → orange → gold → white
    glm::vec3 c = glm::mix(glm::vec3(0.30, 0.03, 0.0), glm::vec3(1.05, 0.16, 0.02), smooth01(0.0, 0.35, t));
    c = glm::mix(c, glm::vec3(1.6, 0.62, 0.10), smooth01(0.30, 0.60, t));
    c = glm::mix(c, glm::vec3(2.0, 1.45, 0.45), smooth01(0.55, 0.82, t));
    c = glm::mix(c, glm::vec3(2.6, 2.35, 2.0), smooth01(0.80, 1.0, t));
    return c;
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 320;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    RigidBody::sleepEpsilon = (real)0.7;

    // --- concrete block, Voronoi-fractured into irregular convex fragments ---
    Vector3 blockHalf(0.85, 0.95, 0.85); Vector3 blockC(0, blockHalf.y + 0.02, 0);
    const real density = 2.4;
    unsigned cs = 4242; auto cr = [&]() { cs = cs * 1103515245u + 12345u; return (float)(((cs >> 16) & 0x7fff) / 32767.0); };
    std::vector<Vector3> seeds;                                            // jittered-grid seeds → varied chunk sizes
    for (int i = 0; i < 78; i++) seeds.push_back(Vector3((cr() * 2 - 1) * blockHalf.x * 0.94, (cr() * 2 - 1) * blockHalf.y * 0.94, (cr() * 2 - 1) * blockHalf.z * 0.94));
    auto cells = voronoiFracture(blockHalf, seeds, 22);

    std::vector<Frag> frags; Destructible dest;
    for (auto& cell : cells) {
        Frag f; f.rad = 0; for (auto& v : cell.verts) f.rad = std::max(f.rad, (float)v.magnitude());
        f.body = std::make_unique<RigidBody>(); f.body->setMass(cell.volume * density);
        Matrix3 I; for (int k = 0; k < 9; k++) I.data[k] = cell.inertiaUnit.data[k] * density; f.body->setInertiaTensor(I);
        f.body->setPosition(blockC + cell.site); f.body->setAcceleration(0, -9.81, 0);
        f.body->setDamping(0.5, 0.5); f.body->calculateDerivedData(); f.body->setAwake(false);   // welded
        f.hull = std::make_unique<CollisionConvex>(); f.hull->body = f.body.get(); f.hull->vertices = cell.verts;
        float t = (cr() - 0.5f) * 0.08f, warm = cr() * 0.04f;                 // concrete tone variation
        f.colour = glm::vec3(0.36f + t + warm, 0.35f + t, 0.33f + t - warm * 0.5f);
        frags.push_back(std::move(f));
    }
    for (auto& f : frags) dest.fragments.push_back(f.body.get());

    // static ground = a big flat convex box (convex–convex contacts, top face at y=0)
    RigidBody gnd; gnd.setInverseMass(0); Matrix3 gi; gi.setInertiaTensorCoeffs(1e12, 1e12, 1e12); gnd.setInertiaTensor(gi);
    gnd.setPosition(Vector3(0, -0.4, 0)); gnd.setDamping(1, 1); gnd.calculateDerivedData();
    CollisionConvex groundHull; groundHull.body = &gnd; groundHull.setBox(Vector3(40, 0.4, 40)); groundHull.calculateInternals();

    World world(8000);
    for (auto& f : frags) world.getRigidBodies().push_back(f.body.get());
    bool frozen = false;                                   // (kept false: contacts always run so pieces rest on locked ones)
    std::vector<char> locked(frags.size(), 0); std::vector<int> pcalm(frags.size(), 0);
    Collider col; col.frags = &frags; col.dest = &dest; col.frozen = &frozen; col.ground = &groundHull;
    col.fr = 0.78; col.re = 0.0;                            // no restitution → debris doesn't keep bouncing
    world.getContactGenerators().push_back(&col);

    // --- blast fire/smoke ---
    std::vector<Puff> puffs;
    Vector3 charge = blockC - Vector3(0, blockHalf.y * 0.35, 0);
    glm::vec3 chargeG((float)charge.x, (float)charge.y, (float)charge.z);
    bool detonated = false; double simT = 0; const double BOOM = 0.35;

    auto detonate = [&]() {
        dest.detonate(charge, 7.5, 0.5, 6.0, 12.0, 11u);
        auto rdir = [&]() { glm::vec3 d(cr() - 0.5f, cr() - 0.5f, cr() - 0.5f); float m = glm::length(d); return m > 1e-4f ? d / m : glm::vec3(0, 1, 0); };
        for (int i = 0; i < 155; i++) {                                       // fireball core + tongues
            Puff q; q.fire = true;
            q.p = chargeG + rdir() * (0.15f + cr() * (float)blockHalf.x * 1.1f);
            q.v = rdir() * (3.0f + cr() * 8.0f) + glm::vec3(0, 2.2f + cr() * 3.2f, 0);   // spread into tongues, rise
            q.life = q.life0 = 0.30f + cr() * 0.6f; q.size = 0.16f + cr() * 0.42f; q.grow = 1.5f + cr() * 1.3f; q.buoy = 4.4f + cr() * 2.8f;
            puffs.push_back(q);
        }
        for (int i = 0; i < 150; i++) {                                       // rising dust/smoke plume (wispier)
            Puff q; q.fire = false;
            q.p = chargeG + rdir() * (0.1f + cr() * (float)blockHalf.x);
            q.v = rdir() * (1.0f + cr() * 2.6f) + glm::vec3(0, 2.2f + cr() * 1.8f, 0);
            q.life = q.life0 = 1.3f + cr() * 2.2f; q.size = 0.18f + cr() * 0.34f; q.grow = 0.8f + cr() * 0.6f; q.buoy = 2.0f + cr() * 1.3f;
            puffs.push_back(q);
        }
    };
    auto stepPhysics = [&]() {
        for (int k = 0; k < 2; k++) { world.startFrame(); world.runPhysics(1.0 / 240); }    // 2 substeps → stabler resting
        simT += 1.0 / 120;
        if (!detonated && simT >= BOOM) { detonate(); detonated = true; }
        else if (detonated) {
            // Single-point convex resting contacts micro-jitter forever, so a pile never
            // fully stops. Two-part cure: (1) bleed the velocity of *slow* pieces each step
            // (fast, still-flying pieces are untouched); (2) once a piece has been slow for a
            // moment, LOCK it into immovable geometry so it can't jitter, while pieces still
            // landing collide and rest on it. The rubble goes fully still, piece by piece.
            double td = simT - BOOM;
            for (size_t i = 0; i < frags.size(); i++) {
                RigidBody* b = frags[i].body.get();
                if (locked[i]) { b->setVelocity(Vector3()); b->setRotation(Vector3()); continue; }
                real v = b->getVelocity().magnitude(), w = b->getRotation().magnitude();
                if (v < 0.6) { b->setVelocity(b->getVelocity() * 0.82); b->setRotation(b->getRotation() * 0.70); }  // angular jitter bled harder
                bool settled = (v < 0.35 && w < 0.9);
                if (settled) pcalm[i]++; else pcalm[i] = 0;
                if (pcalm[i] > 14 || td > 3.4) {            // settled (or past the deadline) → pin it in place
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

    gfx::App app(W, H, "concrete explosion", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh plane = gfx::makePlane(60, 1);
    // build each fragment's flat-shaded convex-polyhedron mesh (needs the GL context)
    for (size_t i = 0; i < frags.size(); i++) {
        auto& cell = cells[i]; std::vector<float> vv; std::vector<unsigned> ii;
        for (size_t t = 0; t < cell.tris.size(); t++) { Vector3 n = cell.triN[t];
            for (int k = 0; k < 3; k++) { Vector3 p = cell.verts[cell.tris[t][k]];
                vv.insert(vv.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z}); ii.push_back((unsigned)ii.size()); } }
        frags[i].mesh.upload(vv, ii);
    }
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 0.9, 0); cam.dist = 7.4f; cam.yaw = 0.7f; cam.pitch = 0.17f;
    glm::mat4 proj = glm::perspective(glm::radians(43.0f), (float)W / H, 0.1f, 100.0f);
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
  float n = vnoise(vUV*2.6+vSeed*4.1)*0.65 + vnoise(vUV*6.0+vSeed*2.3)*0.35;   // 2-octave turbulence
  float soft = smoothstep(1.0,0.12,rr) * (0.12+1.15*n);                        // ragged, wispy tongues
  soft = clamp(soft,0.0,1.0);
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

    auto emitBatch = [&](std::vector<float>& v, const Puff& q, glm::vec3 R, glm::vec3 U) {
        float t = q.life / q.life0; glm::vec4 col;
        if (q.fire) { float temp = glm::clamp(t * 1.15f, 0.0f, 1.0f); glm::vec3 c = fireColor(temp) * (0.35f + 0.85f * temp);
            col = glm::vec4(c, glm::clamp(t * 1.3f, 0.0f, 1.0f)); }
        else { float g = 0.16f + 0.10f * t; float shade = smooth01(0.0f, 0.4f, t); col = glm::vec4(g, g * 0.96f, g * 0.90f, 0.34f * shade); }   // light concrete dust
        const float uv[6][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, -1}, {1, 1}, {-1, 1}};
        for (auto& c : uv) v.insert(v.end(), {q.p.x, q.p.y, q.p.z, c[0], c[1], q.size, col.r, col.g, col.b, col.a});
    };
    auto drawPuffs = [&](glm::vec3 eye, glm::vec3 R, glm::vec3 U) {
        // smoke back-to-front, fire additive on top
        std::vector<float> vs, vf;
        glm::vec3 fwd = glm::normalize(cam.target - eye);
        std::vector<const Puff*> sm; for (auto& q : puffs) if (!q.fire) sm.push_back(&q);
        std::sort(sm.begin(), sm.end(), [&](const Puff* a, const Puff* b) { return glm::dot(a->p - eye, fwd) > glm::dot(b->p - eye, fwd); });
        for (auto* q : sm) emitBatch(vs, *q, R, U);
        for (auto& q : puffs) if (q.fire) emitBatch(vf, q, R, U);
        glUseProgram(pBB);
        glUniformMatrix4fv(glGetUniformLocation(pBB, "uV"), 1, GL_FALSE, &r.curView[0][0]);
        glUniformMatrix4fv(glGetUniformLocation(pBB, "uP"), 1, GL_FALSE, &r.curProj[0][0]);
        glUniform3fv(glGetUniformLocation(pBB, "uRight"), 1, &R[0]); glUniform3fv(glGetUniformLocation(pBB, "uUp"), 1, &U[0]);
        glBindVertexArray(bbVAO); glEnable(GL_BLEND); glDepthMask(GL_FALSE);
        if (!vs.empty()) { glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); glUniform1i(glGetUniformLocation(pBB, "uAdd"), 0);
            glBindBuffer(GL_ARRAY_BUFFER, bbVBO); glBufferData(GL_ARRAY_BUFFER, vs.size() * sizeof(float), vs.data(), GL_DYNAMIC_DRAW);
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(vs.size() / 10)); }
        if (!vf.empty()) { glBlendFunc(GL_ONE, GL_ONE); glUniform1i(glGetUniformLocation(pBB, "uAdd"), 1);
            glBindBuffer(GL_ARRAY_BUFFER, bbVBO); glBufferData(GL_ARRAY_BUFFER, vf.size() * sizeof(float), vf.data(), GL_DYNAMIC_DRAW);
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(vf.size() / 10)); }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
    };

    auto renderFrame = [&]() {
        glm::vec3 eye = cam.eye(); glm::vec3 fwd = glm::normalize(cam.target - eye);
        glm::vec3 R = glm::normalize(glm::cross(fwd, glm::vec3(0, 1, 0))), U = glm::cross(R, fwd);
        r.setLightForScene(glm::vec3(0, 1.0, 0), 7.5f);
        r.beginShadow();
        for (auto& f : frags) r.shadowDraw(f.mesh, glmFromPhys(f.body->getTransform()));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        glClearColor(0.20f, 0.21f, 0.25f, 1.0f); glClear(GL_COLOR_BUFFER_BIT);    // moody demolition sky
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.34, 0.33, 0.31), 0.96f, 0.0f);
        for (auto& f : frags) r.drawPBR(r.pSolid, f.mesh, glmFromPhys(f.body->getTransform()), f.colour, 0.82f, 0.0f);
        drawPuffs(eye, R, U);                                                     // fire + smoke into the HDR buffer
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
