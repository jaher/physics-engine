// Contact-solver upgrade (persistent warm-started manifolds + friction). Tall
// towers of crates stand rock-stable on the improved solver, then a heavy ball is
// launched through them and they topple and scatter, resolved as ordinary contacts.
//   ./stack3d [--shot out.png [f]] | --video f [n]
#include "phys/phys.h"
#include "common/gfx.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <memory>
using namespace phys;

static glm::mat4 glmFromPhys(const Matrix4& t) {
    return glm::mat4(t.data[0], t.data[4], t.data[8], 0, t.data[1], t.data[5], t.data[9], 0,
                     t.data[2], t.data[6], t.data[10], 0, t.data[3], t.data[7], t.data[11], 1);
}
struct Obj { std::unique_ptr<RigidBody> body; bool box; Vector3 half; real radius; glm::vec3 col;
             std::unique_ptr<CollisionBox> cbox; std::unique_ptr<CollisionSphere> csph; };

struct Coll : ContactGenerator {
    std::vector<Obj>* objs; CollisionPlane floor; real fr, re;
    unsigned addContact(Contact* c, unsigned limit) const override {
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
        for (auto& o : *objs) { if (o.box) o.cbox->calculateInternals(); else o.csph->calculateInternals(); }
        for (auto& o : *objs) { if (!d.hasMoreContacts()) break;
            if (o.box) CollisionDetector::boxAndHalfSpace(*o.cbox, floor, &d);
            else CollisionDetector::sphereAndHalfSpace(*o.csph, floor, &d); }
        auto& O = *objs;
        for (size_t i = 0; i < O.size(); i++) for (size_t j = i + 1; j < O.size(); j++) {
            if (!d.hasMoreContacts()) break;
            if ((O[i].body->getPosition() - O[j].body->getPosition()).magnitude() > 1.4) continue;
            if (O[i].box && O[j].box) CollisionDetector::boxAndBox(*O[i].cbox, *O[j].cbox, &d);
            else if (O[i].box) CollisionDetector::boxAndSphere(*O[i].cbox, *O[j].csph, &d);
            else if (O[j].box) CollisionDetector::boxAndSphere(*O[j].cbox, *O[i].csph, &d);
            else CollisionDetector::sphereAndSphere(*O[i].csph, *O[j].csph, &d); }
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

    World world(6000);
    std::vector<Obj> objs;
    unsigned sd = 9u; auto rnd = [&]() { sd = sd * 1103515245u + 12345u; return (float)(((sd >> 16) & 0x7fff) / 32767.0f); };
    auto addBox = [&](Vector3 pos, Vector3 half, glm::vec3 col) {
        Obj o; o.box = true; o.half = half; o.col = col;
        o.body = std::make_unique<RigidBody>(); real m = half.x * half.y * half.z * 40; o.body->setMass(m);
        Matrix3 I; I.setBlockInertiaTensor(half, m); o.body->setInertiaTensor(I);
        o.body->setPosition(pos); o.body->setAcceleration(0, -9.81, 0); o.body->setDamping(0.8, 0.8); o.body->calculateDerivedData();
        o.cbox = std::make_unique<CollisionBox>(); o.cbox->halfSize = half;
        objs.push_back(std::move(o)); objs.back().cbox->body = objs.back().body.get();
        world.getRigidBodies().push_back(objs.back().body.get());
    };
    auto addSphere = [&](Vector3 pos, Vector3 vel, real rad, glm::vec3 col) {
        Obj o; o.box = false; o.radius = rad; o.col = col;
        o.body = std::make_unique<RigidBody>(); real m = rad * rad * rad * 300; o.body->setMass(m);
        Matrix3 I; real cc = 0.4 * m * rad * rad; I.setInertiaTensorCoeffs(cc, cc, cc); o.body->setInertiaTensor(I);
        o.body->setPosition(pos); o.body->setVelocity(vel); o.body->setAcceleration(0, -9.81, 0); o.body->setDamping(0.999, 0.9); o.body->calculateDerivedData();
        o.csph = std::make_unique<CollisionSphere>(); o.csph->radius = rad;
        objs.push_back(std::move(o)); objs.back().csph->body = objs.back().body.get();
        world.getRigidBodies().push_back(objs.back().body.get());
    };
    // a broad brick pyramid — stable on the warm-started manifold solver
    real hs = 0.2, gap = 2 * hs + 0.0004; const int BASE = 6;
    glm::vec3 pal[4] = {{0.80, 0.35, 0.25}, {0.30, 0.55, 0.75}, {0.85, 0.7, 0.25}, {0.4, 0.65, 0.45}};
    for (int row = 0; row < BASE; row++) for (int c = 0; c <= (BASE - 1) - row; c++)
        addBox(Vector3((c - ((BASE - 1) - row) * 0.5) * gap, hs + row * gap, 0), Vector3(hs, hs, hs), pal[(row + c) % 4]);

    Coll coll; coll.objs = &objs; coll.floor.direction = Vector3(0, 1, 0); coll.floor.offset = 0; coll.fr = 0.7; coll.re = 0.0;
    world.getContactGenerators().push_back(&coll);

    gfx::App app(W, H, "stable stacking", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.1, 0.4); cam.dist = 8.2f; cam.yaw = 0.66f; cam.pitch = 0.16f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    int frame = 0;
    auto step = [&]() {
        if (frame == 150) addSphere(Vector3(-5.5, 1.5, 0), Vector3(11, 0.5, 0), 0.38, glm::vec3(0.2, 0.2, 0.22));   // launch the ball
        for (int s = 0; s < 4; s++) { world.startFrame(); world.runPhysics(1.0 / 240); }
        frame++;
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 1.3, 0.4), 6.5f);
        r.beginShadow();
        for (auto& o : objs) { glm::mat4 m = glm::scale(glmFromPhys(o.body->getTransform()), o.box ? glm::vec3(o.half.x, o.half.y, o.half.z) : glm::vec3(o.radius)); r.shadowDraw(o.box ? box : sphere, m); }
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.22, 0.22, 0.24), 0.95f, 0.0f);
        for (auto& o : objs) { glm::mat4 m = glm::scale(glmFromPhys(o.body->getTransform()), o.box ? glm::vec3(o.half.x, o.half.y, o.half.z) : glm::vec3(o.radius));
            r.drawPBR(r.pSolid, o.box ? box : sphere, m, o.col, o.box ? 0.6f : 0.25f, o.box ? 0.0f : 0.8f); }
        r.endScene();
    };

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
