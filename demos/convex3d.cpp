// Convex & mesh geometry (phys::collide_convex): a pile of arbitrary convex hulls
// (faceted gems and cubes) tumbles and stacks via GJK boolean + EPA penetration,
// while analytic cylinders and cones rest on the ground. The same header also has
// a static triangle-mesh collider (sphere/box vs trimesh) — see tests/test_convex.
//   ./convex3d [--shot out.png [f]] | --video f [n]
#include "phys/phys.h"
#include "phys/collide_convex.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <memory>
#include <array>
using namespace phys;

static glm::mat4 glmFromPhys(const Matrix4& t) {
    return glm::mat4(t.data[0], t.data[4], t.data[8], 0, t.data[1], t.data[5], t.data[9], 0,
                     t.data[2], t.data[6], t.data[10], 0, t.data[3], t.data[7], t.data[11], 1);
}
static gfx::Mesh meshFromPoly(const std::vector<Vector3>& v, const std::vector<std::array<int, 3>>& tri) {
    std::vector<float> vd; std::vector<unsigned> idx; unsigned k = 0;
    for (auto& t : tri) {
        Vector3 a = v[t[0]], b = v[t[1]], c = v[t[2]], n = (b - a) % (c - a); if (n.squareMagnitude() > 1e-12) n.normalise();
        for (Vector3 p : {a, b, c}) { vd.insert(vd.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z}); idx.push_back(k++); }
    }
    gfx::Mesh m; m.upload(vd, idx); return m;
}
static gfx::Mesh makeCyl(int seg, real r, real hh, bool cone) {
    std::vector<float> v; std::vector<unsigned> idx;
    auto push = [&](Vector3 p, Vector3 n) { v.insert(v.end(), {(float)p.x, (float)p.y, (float)p.z, (float)n.x, (float)n.y, (float)n.z}); };
    for (int i = 0; i < seg; i++) {
        real a0 = 2 * real_pi * i / seg, a1 = 2 * real_pi * (i + 1) / seg;
        Vector3 b0(r * std::cos(a0), -hh, r * std::sin(a0)), b1(r * std::cos(a1), -hh, r * std::sin(a1));
        Vector3 t0 = cone ? Vector3(0, hh, 0) : Vector3(r * std::cos(a0), hh, r * std::sin(a0));
        Vector3 t1 = cone ? Vector3(0, hh, 0) : Vector3(r * std::cos(a1), hh, r * std::sin(a1));
        Vector3 nn = ((b1 - b0) % (t0 - b0)); if (nn.squareMagnitude() > 1e-9) nn.normalise();
        unsigned s = v.size() / 6; push(b0, nn); push(b1, nn); push(t0, nn); push(t1, nn);
        idx.insert(idx.end(), {s, s + 2, s + 1, s + 1, s + 2, s + 3});
        unsigned c = v.size() / 6; push(Vector3(0, -hh, 0), Vector3(0, -1, 0)); push(b1, Vector3(0, -1, 0)); push(b0, Vector3(0, -1, 0));
        idx.insert(idx.end(), {c, c + 1, c + 2});
    }
    gfx::Mesh m; m.upload(v, idx); return m;
}

struct Obj { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionConvex> hull;
             std::unique_ptr<CollisionCylinder> cyl; std::unique_ptr<CollisionCone> cone;
             int kind; glm::vec3 col; real rad; };   // kind: 0 gem, 1 cube, 2 cyl, 3 cone

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    // gem = hexagonal bipyramid
    std::vector<Vector3> gemV; std::vector<std::array<int, 3>> gemT;
    for (int i = 0; i < 6; i++) { real a = 2 * real_pi * i / 6; gemV.push_back(Vector3(0.16 * std::cos(a), 0, 0.16 * std::sin(a))); }
    gemV.push_back(Vector3(0, 0.22, 0)); gemV.push_back(Vector3(0, -0.16, 0));    // apexes
    for (int i = 0; i < 6; i++) { gemT.push_back({i, 6, (i + 1) % 6}); gemT.push_back({(i + 1) % 6, 7, i}); }
    std::vector<Vector3> cubeV = {{-0.15, -0.15, -0.15}, {0.15, -0.15, -0.15}, {-0.15, 0.15, -0.15}, {0.15, 0.15, -0.15}, {-0.15, -0.15, 0.15}, {0.15, -0.15, 0.15}, {-0.15, 0.15, 0.15}, {0.15, 0.15, 0.15}};

    World world(4000);
    std::vector<Obj> objs;
    // ground = a big flat convex box (for hulls) + a coincident half-space (for cyl/cone)
    RigidBody gnd; gnd.setInverseMass(0); Matrix3 gi; gi.setInertiaTensorCoeffs(1e12, 1e12, 1e12); gnd.setInertiaTensor(gi);
    gnd.setPosition(Vector3(0, -0.3, 0)); gnd.setDamping(1, 1); gnd.calculateDerivedData();
    CollisionConvex groundHull; groundHull.body = &gnd; groundHull.setBox(Vector3(4, 0.3, 4)); groundHull.calculateInternals();
    CollisionPlane floor; floor.direction = Vector3(0, 1, 0); floor.offset = 0;

    unsigned sd = 3u; auto rnd = [&]() { sd = sd * 1103515245u + 12345u; return (float)(((sd >> 16) & 0x7fff) / 32767.0f); };
    auto quatRand = [&]() { Quaternion q(rnd() - 0.5f, rnd() - 0.5f, rnd() - 0.5f, rnd() - 0.5f); q.normalise(); return q; };
    auto addHull = [&](std::vector<Vector3> verts, Vector3 pos, int kind, glm::vec3 col) {
        Obj o; o.kind = kind; o.col = col; o.rad = 0.2;
        o.body = std::make_unique<RigidBody>(); o.body->setMass(1.2); Matrix3 I; I.setBlockInertiaTensor(Vector3(0.16, 0.18, 0.16), 1.2); o.body->setInertiaTensor(I);
        o.body->setPosition(pos); o.body->setOrientation(quatRand()); o.body->setAcceleration(0, -9.81, 0); o.body->setDamping(0.9, 0.85); o.body->calculateDerivedData();
        o.hull = std::make_unique<CollisionConvex>(); o.hull->vertices = verts;
        objs.push_back(std::move(o)); objs.back().hull->body = objs.back().body.get();
        world.getRigidBodies().push_back(objs.back().body.get());
    };
    auto addCylCone = [&](bool cone, Vector3 pos, glm::vec3 col) {
        Obj o; o.kind = cone ? 3 : 2; o.col = col; o.rad = 0.18;
        o.body = std::make_unique<RigidBody>(); o.body->setMass(1.4); Matrix3 I; I.setBlockInertiaTensor(Vector3(0.18, 0.24, 0.18), 1.4); o.body->setInertiaTensor(I);
        o.body->setPosition(pos); o.body->setAcceleration(0, -9.81, 0); o.body->setDamping(0.9, 0.85); o.body->calculateDerivedData();
        if (cone) { o.cone = std::make_unique<CollisionCone>(); o.cone->radius = 0.18; o.cone->halfHeight = 0.24; }
        else { o.cyl = std::make_unique<CollisionCylinder>(); o.cyl->radius = 0.18; o.cyl->halfHeight = 0.24; }
        objs.push_back(std::move(o));
        if (cone) objs.back().cone->body = objs.back().body.get(); else objs.back().cyl->body = objs.back().body.get();
        world.getRigidBodies().push_back(objs.back().body.get());
    };
    glm::vec3 gemCols[4] = {{0.85, 0.25, 0.35}, {0.3, 0.7, 0.85}, {0.9, 0.7, 0.25}, {0.5, 0.8, 0.4}};
    for (int i = 0; i < 14; i++) {                                    // a tumbling pile of gems & cubes
        bool cube = (i % 3 == 0);
        addHull(cube ? cubeV : gemV, Vector3((rnd() - 0.5) * 1.1, 1.4 + i * 0.34, (rnd() - 0.5) * 1.1), cube ? 1 : 0, gemCols[i % 4]);
    }
    addCylCone(false, Vector3(-1.7, 0.9, 0.6), {0.7, 0.55, 0.85}); addCylCone(false, Vector3(1.7, 0.9, -0.4), {0.55, 0.7, 0.8});
    addCylCone(true, Vector3(-1.6, 1.2, -0.7), {0.85, 0.6, 0.4}); addCylCone(true, Vector3(1.6, 1.2, 0.7), {0.6, 0.75, 0.55});

    struct Coll : ContactGenerator {
        std::vector<Obj>* objs; CollisionConvex* ground; CollisionPlane floor; real fr, re;
        unsigned addContact(Contact* c, unsigned limit) const override {
            CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
            for (auto& o : *objs) { if (o.hull) o.hull->calculateInternals(); if (o.cyl) o.cyl->calculateInternals(); if (o.cone) o.cone->calculateInternals(); }
            auto& O = *objs;
            for (auto& o : O) { if (!d.hasMoreContacts()) break;
                if (o.hull) ConvexCollision::convexAndConvex(*o.hull, *ground, &d);
                else if (o.cyl) ConvexCollision::cylinderAndHalfSpace(*o.cyl, floor, &d);
                else if (o.cone) ConvexCollision::coneAndHalfSpace(*o.cone, floor, &d); }
            for (size_t i = 0; i < O.size(); i++) for (size_t j = i + 1; j < O.size(); j++) {
                if (!d.hasMoreContacts()) break;
                if (!O[i].hull || !O[j].hull) continue;              // convex-convex among hulls
                if ((O[i].body->getPosition() - O[j].body->getPosition()).magnitude() > 0.8) continue;
                ConvexCollision::convexAndConvex(*O[i].hull, *O[j].hull, &d); }
            return d.contactCount;
        }
    } coll; coll.objs = &objs; coll.ground = &groundHull; coll.floor = floor; coll.fr = 0.55; coll.re = 0.1;
    world.getContactGenerators().push_back(&coll);

    gfx::App app(W, H, "convex & mesh", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh gemMesh = meshFromPoly(gemV, gemT), boxMesh = gfx::makeBox(), plane = gfx::makePlane(40, 1);
    gfx::Mesh cylMesh = makeCyl(28, 0.18, 0.24, false), coneMesh = makeCyl(28, 0.18, 0.24, true);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 0.7, 0); cam.dist = 5.6f; cam.yaw = 0.6f; cam.pitch = 0.24f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    auto draw = [&](gfx::Renderer& R, Obj& o, bool shadow) {
        glm::mat4 M = glmFromPhys(o.body->getTransform());
        gfx::Mesh* mesh = o.kind == 0 ? &gemMesh : o.kind == 1 ? &boxMesh : o.kind == 2 ? &cylMesh : &coneMesh;
        glm::mat4 m = o.kind == 1 ? glm::scale(M, glm::vec3(0.15f)) : M;
        if (shadow) R.shadowDraw(*mesh, m); else R.drawPBR(R.pSolid, *mesh, m, o.col, 0.3f, 0.15f);
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 0.8, 0), 5.0f);
        r.beginShadow(); for (auto& o : objs) draw(r, o, true); r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.22, 0.22, 0.25), 0.95f, 0.0f);
        glDisable(GL_CULL_FACE);
        for (auto& o : objs) draw(r, o, false);
        r.endScene();
    };
    auto step = [&]() { for (int s = 0; s < 4; s++) { world.startFrame(); world.runPhysics(1.0 / 240); } };

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
