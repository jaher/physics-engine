// A chain reaction: a bullet punches through a glass window, carries on into the
// room, strikes an explosive barrel and detonates it — and the blast shatters the
// crates, pillars and blocks stacked around it. Three fracture models in one shot:
// the pane breaks with an impact ("bullet-hole") pattern (phys::glassfrac) and its
// shards spall/fall; the barrel and the surrounding objects are Voronoi-fractured
// (phys::voronoi) and blasted apart by the explosion. Everything is convex (GJK/EPA),
// contained by the room's wall half-spaces, and settles to rest. Glass is drawn as
// transparent Fresnel-edged material; the blast is an HDR-bloomed fireball + dust.
//   ./bulletbarrel3d --shot out.png [frames]   |   ./bulletbarrel3d --video f [n]
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
static gfx::Mesh makeCylinder(int seg = 34) {
    std::vector<float> v; std::vector<unsigned> idx;
    auto push = [&](float x, float y, float z, float nx, float ny, float nz) { v.insert(v.end(), {x, y, z, nx, ny, nz}); };
    for (int i = 0; i <= seg; i++) { float a = 2 * M_PI * i / seg, c = std::cos(a), s = std::sin(a);
        push(c, 0.5f, s, c, 0, s); push(c, -0.5f, s, c, 0, s); }
    for (int i = 0; i < seg; i++) { unsigned a = i * 2; idx.insert(idx.end(), {a, a + 1, a + 2, a + 2, a + 1, a + 3}); }
    unsigned base = (unsigned)(v.size() / 6); push(0, 0.5f, 0, 0, 1, 0);
    for (int i = 0; i <= seg; i++) { float a = 2 * M_PI * i / seg; push(std::cos(a), 0.5f, std::sin(a), 0, 1, 0); }
    for (int i = 0; i < seg; i++) idx.insert(idx.end(), {base, base + 1 + i, base + 2 + i});
    unsigned bb = (unsigned)(v.size() / 6); push(0, -0.5f, 0, 0, -1, 0);
    for (int i = 0; i <= seg; i++) { float a = 2 * M_PI * i / seg; push(std::cos(a), -0.5f, std::sin(a), 0, -1, 0); }
    for (int i = 0; i < seg; i++) idx.insert(idx.end(), {bb, bb + 2 + i, bb + 1 + i});
    gfx::Mesh m; m.upload(v, idx); return m;
}

struct Frag { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionConvex> hull; gfx::Mesh mesh; float rad; glm::vec3 colour;
              int kind; bool active = false, locked = false; int pcalm = 0; };   // kind: 0 glass, 1 object, 2 barrel

struct Collider : ContactGenerator {
    std::vector<Frag>* fr; const CollisionPlane* walls; int nWalls; real fri, res;
    unsigned addContact(Contact* c, unsigned limit) const override {
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fri; d.restitution = res;
        for (auto& f : *fr) if (f.active) f.hull->calculateInternals();
        for (auto& f : *fr) if (f.active) for (int w = 0; w < nWalls; w++) { if (!d.hasMoreContacts()) break; ConvexCollision::convexAndHalfSpace(*f.hull, walls[w], &d); }
        for (size_t i = 0; i < fr->size(); i++) { if (!(*fr)[i].active) continue; if (!d.hasMoreContacts()) break;
            Vector3 pi = (*fr)[i].body->getPosition(); float ri = (*fr)[i].rad;
            for (size_t j = i + 1; j < fr->size(); j++) { if (!(*fr)[j].active) continue; if (!d.hasMoreContacts()) break;
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
    RigidBody::sleepEpsilon = (real)0.6;

    // room
    const real Rx = 6.0, Rz = 6.5, Hy = 5.0;
    CollisionPlane walls[6];
    walls[0].direction = Vector3(0, 1, 0);  walls[0].offset = 0;
    walls[1].direction = Vector3(0, -1, 0); walls[1].offset = -Hy;
    walls[2].direction = Vector3(-1, 0, 0); walls[2].offset = -Rx;
    walls[3].direction = Vector3(1, 0, 0);  walls[3].offset = -Rx;
    walls[4].direction = Vector3(0, 0, -1); walls[4].offset = -Rz;
    walls[5].direction = Vector3(0, 0, 1);  walls[5].offset = -Rz;

    std::vector<Frag> frags; std::vector<ConvexCell> cells;

    // glass window (kind 0)
    const glm::vec3 paneC(0, 1.1f, -1.8f); const real PW = 1.4, PH = 1.0, PT = 0.03;
    const real iu = 0.1, iv = -0.2; glm::vec3 impact(paneC.x + (float)iu, paneC.y + (float)iv, paneC.z);
    for (auto& cell : glassFracture(PW, PH, PT, iu, iv, 20, 7, 0.05, 7u)) {
        Frag f; f.kind = 0; f.rad = 0; for (auto& v : cell.verts) f.rad = std::max(f.rad, (float)v.magnitude());
        f.body = std::make_unique<RigidBody>(); f.body->setMass(std::max(1e-3, cell.volume * 2.5));
        Matrix3 I; for (int k = 0; k < 9; k++) I.data[k] = cell.inertiaUnit.data[k] * 2.5; for (int k : {0, 4, 8}) I.data[k] = std::max(I.data[k], (real)1e-6); f.body->setInertiaTensor(I);
        f.body->setPosition(Vector3(paneC.x + cell.site.x, paneC.y + cell.site.y, paneC.z + cell.site.z));
        f.body->setAcceleration(0, -9.81, 0); f.body->setDamping(0.6, 0.4); f.body->calculateDerivedData(); f.body->setAwake(false);
        f.hull = std::make_unique<CollisionConvex>(); f.hull->body = f.body.get(); f.hull->vertices = cell.verts;
        f.colour = glm::vec3(0.55, 0.72, 0.82);
        frags.push_back(std::move(f)); cells.push_back(std::move(cell));
    }

    // Voronoi object / barrel builder
    auto addChunks = [&](Vector3 center, Vector3 half, int nseed, real density, glm::vec3 base, int kind, unsigned seed) {
        unsigned s = seed; auto rr = [&]() { s = s * 1103515245u + 12345u; return (float)(((s >> 16) & 0x7fff) / 32767.0); };
        std::vector<Vector3> seeds;
        for (int i = 0; i < nseed; i++) seeds.push_back(Vector3((rr() * 2 - 1) * half.x * 0.9, (rr() * 2 - 1) * half.y * 0.9, (rr() * 2 - 1) * half.z * 0.9));
        for (auto& cell : voronoiFracture(half, seeds, 20)) {
            Frag f; f.kind = kind; f.rad = 0; for (auto& v : cell.verts) f.rad = std::max(f.rad, (float)v.magnitude());
            f.body = std::make_unique<RigidBody>(); f.body->setMass(cell.volume * density);
            Matrix3 I; for (int k = 0; k < 9; k++) I.data[k] = cell.inertiaUnit.data[k] * density; f.body->setInertiaTensor(I);
            f.body->setPosition(center + cell.site); f.body->setAcceleration(0, -9.81, 0);
            f.body->setDamping(0.5, 0.5); f.body->calculateDerivedData(); f.body->setAwake(false);
            f.hull = std::make_unique<CollisionConvex>(); f.hull->body = f.body.get(); f.hull->vertices = cell.verts;
            float t = (rr() - 0.5f) * 0.08f; f.colour = base + glm::vec3(t, t * 0.9f, t * 0.8f);
            frags.push_back(std::move(f)); cells.push_back(std::move(cell));
        }
    };
    // barrel (kind 2) + surrounding objects (kind 1)
    const glm::vec3 barrelC(0.1f, 0.7f, 1.4f); const float barrelR = 0.42f, barrelHh = 0.7f;
    {   // hollow barrel: its thin wall + end caps dice into small curved metal shards
        unsigned bs = 55u; auto br = [&]() { bs = bs * 1103515245u + 12345u; return (float)(((bs >> 16) & 0x7fff) / 32767.0); };
        for (auto& cell : cylinderShellFracture(barrelR, barrelHh, 0.03, 15, 5, true, 55u)) {
            Frag f; f.kind = 2; f.rad = 0; for (auto& v : cell.verts) f.rad = std::max(f.rad, (float)v.magnitude());
            f.body = std::make_unique<RigidBody>(); f.body->setMass(std::max(1e-3, cell.volume * 5.0));
            Matrix3 I; for (int k = 0; k < 9; k++) I.data[k] = cell.inertiaUnit.data[k] * 5.0; for (int k : {0, 4, 8}) I.data[k] = std::max(I.data[k], (real)1e-6); f.body->setInertiaTensor(I);
            f.body->setPosition(Vector3(barrelC.x + cell.site.x, barrelC.y + cell.site.y, barrelC.z + cell.site.z));
            f.body->setAcceleration(0, -9.81, 0); f.body->setDamping(0.55, 0.5); f.body->calculateDerivedData(); f.body->setAwake(false);
            f.hull = std::make_unique<CollisionConvex>(); f.hull->body = f.body.get(); f.hull->vertices = cell.verts;
            f.colour = (br() < 0.55f) ? glm::vec3(0.50, 0.17, 0.14) + glm::vec3((br() - 0.5f) * 0.06f)    // red paint
                                      : glm::vec3(0.34, 0.34, 0.37) + glm::vec3((br() - 0.5f) * 0.05f);   // bare torn metal
            frags.push_back(std::move(f)); cells.push_back(std::move(cell));
        }
    }
    glm::vec3 grey(0.52, 0.51, 0.49), sand(0.60, 0.49, 0.31), slate(0.42, 0.43, 0.47);
    addChunks(Vector3(-1.9, 0.6, 2.1),  Vector3(0.60, 0.60, 0.60), 26, 2.4, grey, 1, 13u);
    addChunks(Vector3(1.9, 1.15, 2.4),  Vector3(0.36, 1.15, 0.36), 22, 2.2, sand, 1, 21u);
    addChunks(Vector3(-1.3, 0.55, 3.5), Vector3(0.55, 0.55, 0.55), 22, 2.3, grey, 1, 29u);
    addChunks(Vector3(0.5, 0.4, 3.7),   Vector3(1.00, 0.40, 0.65), 28, 2.3, slate, 1, 41u);
    addChunks(Vector3(2.9, 1.05, 3.3),  Vector3(0.34, 1.05, 0.34), 20, 2.2, sand, 1, 47u);

    World world(24000);
    for (auto& f : frags) world.getRigidBodies().push_back(f.body.get());
    Collider col; col.fr = &frags; col.walls = walls; col.nWalls = 6; col.fri = 0.6; col.res = 0.12;
    world.getContactGenerators().push_back(&col);

    // effects + bullet
    std::vector<Puff> puffs; unsigned pcs = 91u; auto cr = [&]() { pcs = pcs * 1103515245u + 12345u; return (float)(((pcs >> 16) & 0x7fff) / 32767.0); };
    glm::vec3 bpos(impact.x, impact.y, -3.0f); glm::vec3 bvel(0, 0, 16.0f); const float bR = 0.05f;
    bool glassBroke = false, boom = false; double simT = 0, boomT = -1;
    const float glassZ = paneC.z, barrelHitZ = barrelC.z - barrelR + 0.14f;   // detonate once slightly inside the drum
    // slow-motion "bullet time" (ramps back to real speed after the blast) + spiral trail
    double timeScale = 0.16; float bulletDist = 0, spinPhase = 0; const float trailR = 0.04f, spinFreq = 17.0f;   // helix ≈ bullet radius
    std::vector<glm::vec4> trail;   // xyz = corkscrew point, w = brightness

    auto breakGlass = [&]() {
        unsigned st = 3u; auto rr = [&]() { st = st * 1103515245u + 12345u; return (real)(((st >> 16) & 0x7fff) / 32767.0); };
        for (auto& f : frags) if (f.kind == 0) {
            Vector3 p = f.body->getPosition(); real du = p.x - impact.x, dv = p.y - impact.y, dist = std::sqrt(du * du + dv * dv);
            Vector3 rad = dist > 1e-4 ? Vector3(du, dv, 0) * (1.0 / dist) : Vector3(0, 1, 0);
            Vector3 v = Vector3(0, 0, 1) * (2.6 / (dist + 0.18)) + rad * (1.0 / (dist + 0.25)) + Vector3(rr() - 0.5, rr() - 0.5, rr() - 0.5) * 0.6;
            f.active = true; f.body->setAwake(true); f.body->setVelocity(v);
            f.body->setRotation(Vector3(rr() - 0.5, rr() - 0.5, rr() - 0.5) * (4.0 + 8.0 / (dist + 0.3)));
        }
    };
    auto detonate = [&]() {
        unsigned st = 9u; auto rr = [&]() { st = st * 1103515245u + 12345u; return (real)(((st >> 16) & 0x7fff) / 32767.0); };
        Vector3 bc(barrelC.x, barrelC.y, barrelC.z);
        for (auto& f : frags) {
            Vector3 dr = f.body->getPosition() - bc; real dist = dr.magnitude();
            Vector3 dir = dist > 1e-6 ? dr * (1.0 / dist) : Vector3(0, 1, 0);
            real speed = (f.kind == 2 ? 13.0 : 11.0) / (dist + 0.5);
            Vector3 blast = dir * speed + Vector3(0, 0.45 * speed + 1.5, 0);
            if (f.kind == 0) { f.body->setVelocity(f.body->getVelocity() + blast * 0.5); }   // push falling glass
            else { f.active = true; f.body->setAwake(true); f.body->setVelocity(blast + Vector3(rr() - 0.5, rr() - 0.5, rr() - 0.5) * 2.0);
                   f.body->setRotation(Vector3(rr() - 0.5, rr() - 0.5, rr() - 0.5) * 16.0); }
        }
        auto rdir = [&]() { glm::vec3 d(cr() - 0.5f, cr() - 0.5f, cr() - 0.5f); float m = glm::length(d); return m > 1e-4f ? d / m : glm::vec3(0, 1, 0); };
        for (int i = 0; i < 240; i++) { Puff q; q.fire = true; q.p = barrelC + rdir() * (0.2f + cr() * 1.3f);
            q.v = rdir() * (4.0f + cr() * 11.0f) + glm::vec3(0, 2.5f + cr() * 3.5f, 0);
            q.life = q.life0 = 0.32f + cr() * 0.7f; q.size = 0.22f + cr() * 0.55f; q.grow = 1.8f + cr() * 1.6f; q.buoy = 4.6f + cr() * 3.0f; puffs.push_back(q); }
        for (int i = 0; i < 200; i++) { Puff q; q.fire = false; q.p = barrelC + rdir() * (0.15f + cr() * 1.2f);
            q.v = rdir() * (1.4f + cr() * 3.5f) + glm::vec3(0, 2.4f + cr() * 2.0f, 0);
            q.life = q.life0 = 1.6f + cr() * 2.4f; q.size = 0.24f + cr() * 0.45f; q.grow = 0.9f + cr() * 0.7f; q.buoy = 2.0f + cr() * 1.4f; puffs.push_back(q); }
    };
    auto stepPhysics = [&]() {
        if (boom) timeScale = std::min(1.0, timeScale + 0.035);              // ramp back to real speed after the blast
        double simDt = timeScale * (1.0 / 120.0); simT += simDt;
        for (int k = 0; k < 2; k++) { world.startFrame(); world.runPhysics(simDt * 0.5); }
        if (!boom) {
            glm::vec3 stepv = bvel * (float)simDt; bpos += stepv; bulletDist += glm::length(stepv);
            spinPhase = bulletDist * spinFreq;                              // corkscrew: helix around the flight path
            glm::vec3 h = glm::vec3(0, 1, 0) * std::cos(spinPhase) + glm::vec3(1, 0, 0) * std::sin(spinPhase);
            trail.push_back(glm::vec4(bpos + h * trailR, 1.0f));            // double-helix + a faint core
            trail.push_back(glm::vec4(bpos - h * trailR, 1.0f));
            trail.push_back(glm::vec4(bpos, 0.7f));
            if (!glassBroke && bpos.z >= glassZ) { breakGlass(); glassBroke = true; bvel *= 0.85f; }
            if (bpos.z >= barrelHitZ) { detonate(); boom = true; boomT = simT; }
        }
        float tfade = boom ? 0.90f : 0.996f;                                // trail lingers in flight, burns off in the blast
        for (auto& t : trail) t.w *= tfade;
        trail.erase(std::remove_if(trail.begin(), trail.end(), [](const glm::vec4& t) { return t.w < 0.04f; }), trail.end());
        if (trail.size() > 1500) trail.erase(trail.begin(), trail.begin() + (trail.size() - 1500));
        for (auto& f : frags) {                                                 // settle: bleed slow, lock at rest
            if (!f.active || f.locked) { if (f.locked) { f.body->setVelocity(Vector3()); f.body->setRotation(Vector3()); } continue; }
            real v = f.body->getVelocity().magnitude(), w = f.body->getRotation().magnitude();
            if (v < 0.55) { f.body->setVelocity(f.body->getVelocity() * 0.83); f.body->setRotation(f.body->getRotation() * 0.72); }
            if (v < 0.32 && w < 0.85) f.pcalm++; else f.pcalm = 0;
            if (f.pcalm > 14 || (boom && simT - boomT > 3.5)) { f.body->setVelocity(Vector3()); f.body->setRotation(Vector3()); f.body->setAcceleration(0, 0, 0);
                f.body->setInverseMass(0); Matrix3 bi; bi.setInertiaTensorCoeffs(1e12, 1e12, 1e12); f.body->setInertiaTensor(bi); f.body->calculateDerivedData(); f.body->setAwake(false); f.locked = true; }
        }
        float fdt = (float)simDt;
        for (auto& q : puffs) { q.v.y += q.buoy * fdt; q.v *= (q.fire ? 0.90f : 0.95f); q.p += q.v * fdt; q.size += q.grow * fdt; q.life -= fdt; }
        puffs.erase(std::remove_if(puffs.begin(), puffs.end(), [](const Puff& q) { return q.life <= 0; }), puffs.end());
    };

    gfx::App app(W, H, "bullet → barrel → blast", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(24, 1), cyl = makeCylinder();
    for (size_t i = 0; i < frags.size(); i++) { auto& cell = cells[i]; std::vector<float> vv; std::vector<unsigned> ii;
        for (size_t t = 0; t < cell.tris.size(); t++) { Vector3 n = cell.triN[t];
            for (int k = 0; k < 3; k++) { Vector3 p = cell.verts[cell.tris[t][k]]; vv.insert(vv.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z}); ii.push_back((unsigned)ii.size()); } }
        frags[i].mesh.upload(vv, ii); }
    gfx::Mesh pane; { std::vector<float> v = {-PW, -PH, 0, 0, 0, -1, PW, -PH, 0, 0, 0, -1, PW, PH, 0, 0, 0, -1, -PW, PH, 0, 0, 0, -1};
        std::vector<unsigned> idx = {0, 1, 2, 0, 2, 3}; pane.upload(v, idx); }

    gfx::OrbitCamera cam; cam.target = glm::vec3(0.0, 1.0, 0.5); cam.dist = 7.6f; cam.yaw = 2.32f; cam.pitch = 0.22f;
    glm::mat4 proj = glm::perspective(glm::radians(43.0f), (float)W / H, 0.05f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    GLuint pBB = gfx::program(R"(#version 330 core
layout(location=0) in vec3 aC; layout(location=1) in vec2 aUV; layout(location=2) in float aS; layout(location=3) in vec4 aCol;
uniform mat4 uV,uP; uniform vec3 uRight,uUp; out vec2 vUV; out vec4 vCol; out vec2 vSeed;
void main(){ vUV=aUV; vCol=aCol; vSeed=aC.xz; vec3 wp=aC+(aUV.x*uRight+aUV.y*uUp)*aS; gl_Position=uP*uV*vec4(wp,1.0); })",
        R"(#version 330 core
in vec2 vUV; in vec4 vCol; in vec2 vSeed; out vec4 frag; uniform int uAdd;
float hash(vec2 p){ return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453); }
float vnoise(vec2 p){ vec2 i=floor(p),f=fract(p); f=f*f*(3.0-2.0*f);
  float a=hash(i),b=hash(i+vec2(1,0)),c=hash(i+vec2(0,1)),d=hash(i+vec2(1,1)); return mix(mix(a,b,f.x),mix(c,d,f.x),f.y); }
void main(){ float rr=length(vUV); if(rr>1.0) discard;
  float n=vnoise(vUV*2.6+vSeed*4.1)*0.65+vnoise(vUV*6.0+vSeed*2.3)*0.35;
  float soft=smoothstep(1.0,0.12,rr)*(0.12+1.15*n); soft=clamp(soft,0.0,1.0); if(soft<0.06) discard;
  if(uAdd==1) frag=vec4(vCol.rgb*soft*vCol.a,1.0); else frag=vec4(vCol.rgb,soft*vCol.a); })");
    GLuint bbVAO, bbVBO; glGenVertexArrays(1, &bbVAO); glGenBuffers(1, &bbVBO); glBindVertexArray(bbVAO); glBindBuffer(GL_ARRAY_BUFFER, bbVBO);
    glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1); glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2); glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)(5 * sizeof(float)));
    glEnableVertexAttribArray(3); glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 10 * sizeof(float), (void*)(6 * sizeof(float))); glBindVertexArray(0);
    GLuint pGlass = gfx::program(R"(#version 330 core
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aN; uniform mat4 uM,uV,uP; out vec3 vW; out vec3 vN;
void main(){ vec4 w=uM*vec4(aPos,1.0); vW=w.xyz; vN=mat3(uM)*aN; gl_Position=uP*uV*w; })",
        R"(#version 330 core
in vec3 vW; in vec3 vN; out vec4 frag; uniform vec3 uCam,uLight;
void main(){ vec3 N=normalize(vN); vec3 V=normalize(uCam-vW); if(dot(N,V)<0.0) N=-N; vec3 L=normalize(uLight);
  float fres=pow(1.0-abs(dot(N,V)),3.0); float spec=pow(max(dot(reflect(-L,N),V),0.0),80.0); vec3 tint=vec3(0.55,0.73,0.82);
  vec3 col=tint*(0.35+0.35*max(dot(N,L),0.0))+vec3(0.8,0.9,1.0)*fres*0.9+vec3(spec)*1.2; float a=clamp(0.12+0.82*fres+spec,0.0,0.92); frag=vec4(col,a); })");
    auto setM = [](GLuint p, const char* n, const glm::mat4& m) { glUniformMatrix4fv(glGetUniformLocation(p, n), 1, GL_FALSE, &m[0][0]); };
    auto setV = [](GLuint p, const char* n, const glm::vec3& v) { glUniform3fv(glGetUniformLocation(p, n), 1, &v[0]); };

    auto drawPuffs = [&](glm::vec3 eye, glm::vec3 R, glm::vec3 U) {
        std::vector<float> vs, vf; glm::vec3 fwd = glm::normalize(cam.target - eye);
        auto emit = [&](std::vector<float>& out, const Puff& q) { float t = q.life / q.life0; glm::vec4 col;
            if (q.fire) { float tp = glm::clamp(t * 1.15f, 0.0f, 1.0f); glm::vec3 c = fireColor(tp) * (0.35f + 0.85f * tp); col = glm::vec4(c, glm::clamp(t * 1.3f, 0.0f, 1.0f)); }
            else { float g = 0.16f + 0.10f * t; col = glm::vec4(g, g * 0.96f, g * 0.9f, 0.34f * smooth01(0.0f, 0.4f, t)); }
            const float uv[6][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, -1}, {1, 1}, {-1, 1}};
            for (auto& c : uv) out.insert(out.end(), {q.p.x, q.p.y, q.p.z, c[0], c[1], q.size, col.r, col.g, col.b, col.a}); };
        std::vector<const Puff*> sm; for (auto& q : puffs) if (!q.fire) sm.push_back(&q);
        std::sort(sm.begin(), sm.end(), [&](const Puff* a, const Puff* b) { return glm::dot(a->p - eye, fwd) > glm::dot(b->p - eye, fwd); });
        for (auto* q : sm) emit(vs, *q); for (auto& q : puffs) if (q.fire) emit(vf, q);
        glUseProgram(pBB); setM(pBB, "uV", r.curView); setM(pBB, "uP", r.curProj); setV(pBB, "uRight", R); setV(pBB, "uUp", U);
        glBindVertexArray(bbVAO); glEnable(GL_BLEND); glDepthMask(GL_FALSE);
        if (!vs.empty()) { glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); glUniform1i(glGetUniformLocation(pBB, "uAdd"), 0);
            glBindBuffer(GL_ARRAY_BUFFER, bbVBO); glBufferData(GL_ARRAY_BUFFER, vs.size() * sizeof(float), vs.data(), GL_DYNAMIC_DRAW); glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(vs.size() / 10)); }
        if (!vf.empty()) { glBlendFunc(GL_ONE, GL_ONE); glUniform1i(glGetUniformLocation(pBB, "uAdd"), 1);
            glBindBuffer(GL_ARRAY_BUFFER, bbVBO); glBufferData(GL_ARRAY_BUFFER, vf.size() * sizeof(float), vf.data(), GL_DYNAMIC_DRAW); glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(vf.size() / 10)); }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
    };

    struct WallDraw { glm::vec3 c, s, nIn, col; };
    std::vector<WallDraw> wd = {
        {{Rx + 0.1f, Hy * 0.5f, 0}, {0.1f, Hy * 0.5f, Rz}, {-1, 0, 0}, {0.5, 0.5, 0.53}}, {{-Rx - 0.1f, Hy * 0.5f, 0}, {0.1f, Hy * 0.5f, Rz}, {1, 0, 0}, {0.5, 0.5, 0.53}},
        {{0, Hy * 0.5f, Rz + 0.1f}, {Rx, Hy * 0.5f, 0.1f}, {0, 0, -1}, {0.53, 0.52, 0.5}}, {{0, Hy * 0.5f, -Rz - 0.1f}, {Rx, Hy * 0.5f, 0.1f}, {0, 0, 1}, {0.53, 0.52, 0.5}},
        {{0, Hy + 0.1f, 0}, {Rx, 0.1f, Rz}, {0, -1, 0}, {0.42, 0.42, 0.46}}};

    auto renderFrame = [&]() {
        glm::vec3 eye = cam.eye(); glm::vec3 fwd = glm::normalize(cam.target - eye);
        glm::vec3 R = glm::normalize(glm::cross(fwd, glm::vec3(0, 1, 0))), U = glm::cross(R, fwd);
        r.setLightForScene(glm::vec3(0, 1.4, 0.8), 8.0f);
        r.beginShadow();
        for (auto& f : frags) if (f.kind != 0 && (f.active || f.kind == 1)) r.shadowDraw(f.mesh, glmFromPhys(f.body->getTransform()));
        if (!boom) r.shadowDraw(cyl, glm::scale(glm::translate(glm::mat4(1), barrelC), glm::vec3(barrelR, barrelHh * 2, barrelR)));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        glClearColor(0.16f, 0.17f, 0.20f, 1.0f); glClear(GL_COLOR_BUFFER_BIT);
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.30, 0.29, 0.28), 0.96f, 0.0f);
        for (auto& w : wd) if (glm::dot(eye - w.c, w.nIn) > 0) r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), w.c), w.s), w.col, 0.92f, 0.0f);
        // window frame
        glm::vec3 fc(0.11, 0.09, 0.08); float fb = 0.09f;
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(paneC.x, paneC.y + PH + fb, paneC.z)), glm::vec3(PW + fb, fb, 0.08)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(paneC.x, paneC.y - PH - fb, paneC.z)), glm::vec3(PW + fb, fb, 0.08)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(paneC.x - PW - fb, paneC.y, paneC.z)), glm::vec3(fb, PH, 0.08)), fc, 0.7f, 0.1f);
        r.drawPBR(r.pSolid, box, glm::scale(glm::translate(glm::mat4(1), glm::vec3(paneC.x + PW + fb, paneC.y, paneC.z)), glm::vec3(fb, PH, 0.08)), fc, 0.7f, 0.1f);
        // opaque solids: object fragments (rough), barrel shards (metallic), and the intact barrel before the boom
        for (auto& f : frags) if (f.kind != 0 && (f.kind == 1 || f.active)) {
            float rough = f.kind == 2 ? 0.42f : 0.8f, metal = f.kind == 2 ? 0.65f : 0.0f;
            r.drawPBR(r.pSolid, f.mesh, glmFromPhys(f.body->getTransform()), f.colour, rough, metal);
        }
        if (!boom) {   // an open-top hollow drum: red wall, band, dark recessed interior
            r.drawPBR(r.pSolid, cyl, glm::scale(glm::translate(glm::mat4(1), barrelC), glm::vec3(barrelR, barrelHh * 2, barrelR)), glm::vec3(0.55, 0.15, 0.11), 0.45f, 0.4f);
            r.drawPBR(r.pSolid, cyl, glm::scale(glm::translate(glm::mat4(1), barrelC + glm::vec3(0, barrelHh * 0.4f, 0)), glm::vec3(barrelR * 1.02f, 0.05f, barrelR * 1.02f)), glm::vec3(0.10, 0.09, 0.09), 0.6f, 0.3f);
            r.drawPBR(r.pSolid, cyl, glm::scale(glm::translate(glm::mat4(1), barrelC + glm::vec3(0, barrelHh - 0.03f, 0)), glm::vec3(barrelR * 0.90f, 0.05f, barrelR * 0.90f)), glm::vec3(0.08, 0.07, 0.07), 0.85f, 0.0f); }   // hollow opening
        r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), bpos), glm::vec3(bR)), glm::vec3(0.05, 0.05, 0.06), 0.3f, 0.8f);   // bullet

        drawPuffs(eye, R, U);
        // spiral bullet trail (additive glowing corkscrew)
        if (!trail.empty()) {
            std::vector<float> tv; tv.reserve(trail.size() * 60);
            const float uv[6][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, -1}, {1, 1}, {-1, 1}};
            for (auto& t : trail) { glm::vec3 col = glm::vec3(1.0f, 0.82f, 0.5f) * (t.w * 0.85f);
                for (auto& c : uv) tv.insert(tv.end(), {t.x, t.y, t.z, c[0], c[1], 0.016f, col.r, col.g, col.b, 1.0f}); }
            glUseProgram(pBB); setM(pBB, "uV", r.curView); setM(pBB, "uP", r.curProj); setV(pBB, "uRight", R); setV(pBB, "uUp", U);
            glUniform1i(glGetUniformLocation(pBB, "uAdd"), 1);
            glBindVertexArray(bbVAO); glEnable(GL_BLEND); glBlendFunc(GL_ONE, GL_ONE); glDepthMask(GL_FALSE);
            glBindBuffer(GL_ARRAY_BUFFER, bbVBO); glBufferData(GL_ARRAY_BUFFER, tv.size() * sizeof(float), tv.data(), GL_DYNAMIC_DRAW);
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(tv.size() / 10));
            glDepthMask(GL_TRUE); glDisable(GL_BLEND);
        }
        // glass pass
        glUseProgram(pGlass); setM(pGlass, "uV", r.curView); setM(pGlass, "uP", r.curProj); setV(pGlass, "uCam", eye); setV(pGlass, "uLight", glm::normalize(glm::vec3(-0.4, 0.7, 0.5)));
        glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); glDepthMask(GL_FALSE);
        if (!glassBroke) { setM(pGlass, "uM", glm::translate(glm::mat4(1), glm::vec3(paneC))); pane.draw(); }
        else { std::vector<int> gi; for (int i = 0; i < (int)frags.size(); i++) if (frags[i].kind == 0) gi.push_back(i);
            std::sort(gi.begin(), gi.end(), [&](int a, int b) { return glm::length(glm::vec3(frags[a].body->getPosition().x, frags[a].body->getPosition().y, frags[a].body->getPosition().z) - eye) > glm::length(glm::vec3(frags[b].body->getPosition().x, frags[b].body->getPosition().y, frags[b].body->getPosition().z) - eye); });
            for (int i : gi) { setM(pGlass, "uM", glmFromPhys(frags[i].body->getTransform())); frags[i].mesh.draw(); } }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);
        r.endScene();
    };

    if (video) { for (int f = 0; f < frames; f++) { cam.yaw += 0.45f / frames; renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p); for (int s = 0; s < 2; s++) stepPhysics(); }
        std::printf("wrote %d frames (%zu frags)\n", frames, frags.size()); return 0; }
    if (headless) { for (int i = 0; i < frames * 2; i++) stepPhysics(); renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) { if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my); bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down; stepPhysics(); renderFrame(); r.present(); app.poll(); }
    return 0;
}
