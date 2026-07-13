// Conveyor + fields (phys::constraint2). Crates spawn at the left, a ConveyorSurface
// drags them along a moving belt, a crosswind (WindField) nudges them, and they
// tumble off the end into a pile — all on the engine's rigid-body collision resolver.
//   ./constraints3d [--shot out.png [f]] | --video f [n]
#include "phys/phys.h"
#include "phys/constraint2.h"
#include "common/gfx.h"
#include <vector>
#include <memory>
#include <cstring>
#include <cmath>
using namespace phys;

static glm::mat4 glmFromPhys(const Matrix4& t) {
    return glm::mat4(t.data[0], t.data[4], t.data[8], 0, t.data[1], t.data[5], t.data[9], 0,
                     t.data[2], t.data[6], t.data[10], 0, t.data[3], t.data[7], t.data[11], 1);
}
struct Crate { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionBox> cbox; Vector3 half; glm::vec3 col; };

struct Scene : ContactGenerator {
    std::vector<Crate>* crates; CollisionBox* belt; CollisionPlane floor; real fr, re;
    unsigned addContact(Contact* c, unsigned limit) const override {
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
        for (auto& o : *crates) o.cbox->calculateInternals();
        belt->calculateInternals();
        for (auto& o : *crates) { if (!d.hasMoreContacts()) break;
            CollisionDetector::boxAndHalfSpace(*o.cbox, floor, &d);
            CollisionDetector::boxAndBox(*belt, *o.cbox, &d); }
        auto& C = *crates;
        for (size_t i = 0; i < C.size(); i++) for (size_t j = i + 1; j < C.size(); j++) {
            if (!d.hasMoreContacts()) break;
            if ((C[i].body->getPosition() - C[j].body->getPosition()).magnitude() > 1.0) continue;
            CollisionDetector::boxAndBox(*C[i].cbox, *C[j].cbox, &d); }
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
    const real beltY = 1.25, beltHalfX = 2.1, beltHalfZ = 0.55, beltTop = beltY + 0.06;

    World world(3000);
    std::vector<Crate> crates;
    RigidBody beltBody; beltBody.setInverseMass(0); Matrix3 bi; bi.setInertiaTensorCoeffs(1e12, 1e12, 1e12); beltBody.setInertiaTensor(bi);
    beltBody.setPosition(Vector3(0, beltY, 0)); beltBody.setDamping(1, 1); beltBody.calculateDerivedData();
    CollisionBox beltBox; beltBox.body = &beltBody; beltBox.halfSize = Vector3(beltHalfX, 0.06, beltHalfZ); beltBox.calculateInternals();

    Scene scene; scene.crates = &crates; scene.belt = &beltBox;
    scene.floor.direction = Vector3(0, 1, 0); scene.floor.offset = 0; scene.fr = 0.6; scene.re = 0.1;
    world.getContactGenerators().push_back(&scene);

    ConveyorSurface belt; belt.normal = Vector3(0, 1, 0); belt.offset = beltTop; belt.beltVelocity = Vector3(1.4, 0, 0); belt.maxHeight = 0.35; belt.mu = 1.2;
    WindField wind; wind.wind = Vector3(0, 0, 0); wind.gustAxis = Vector3(0, 0, 1); wind.gustAmp = 1.6; wind.gustFreq = 1.1; wind.linearDrag = 0.35; wind.quadraticDrag = 0.05;

    unsigned sd = 5u; auto rnd = [&]() { sd = sd * 1103515245u + 12345u; return (float)(((sd >> 16) & 0x7fff) / 32767.0f); };
    auto spawn = [&]() {
        Crate o; real h = 0.13f + 0.05f * rnd(); o.half = Vector3(h, h, h * (0.8 + 0.4 * rnd())); o.col = glm::vec3(0.35 + 0.5 * rnd(), 0.3 + 0.45 * rnd(), 0.3 + 0.45 * rnd());
        o.body = std::make_unique<RigidBody>(); real m = o.half.x * o.half.y * o.half.z * 60; o.body->setMass(m);
        Matrix3 I; I.setBlockInertiaTensor(o.half, m); o.body->setInertiaTensor(I);
        o.body->setPosition(Vector3(-beltHalfX + 0.3, beltTop + o.half.y + 0.05, (rnd() - 0.5) * 2 * (beltHalfZ - 0.25)));
        o.body->setAcceleration(0, -9.81, 0); o.body->setDamping(0.9, 0.8); o.body->calculateDerivedData();
        o.cbox = std::make_unique<CollisionBox>(); o.cbox->halfSize = o.half;
        crates.push_back(std::move(o)); crates.back().cbox->body = crates.back().body.get();
        world.getRigidBodies().push_back(crates.back().body.get());
    };

    gfx::App app(W, H, "conveyor + fields", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), plane = gfx::makePlane(40, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(-0.2, 0.9, 0); cam.dist = 6.2f; cam.yaw = 0.66f; cam.pitch = 0.24f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });
    glm::mat4 beltM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, beltY, 0)), glm::vec3(beltHalfX, 0.06f, beltHalfZ));
    glm::mat4 legL = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-beltHalfX + 0.2f, beltY * 0.5f, 0)), glm::vec3(0.07f, beltY * 0.5f, beltHalfZ));
    glm::mat4 legR = glm::scale(glm::translate(glm::mat4(1), glm::vec3(beltHalfX - 0.2f, beltY * 0.5f, 0)), glm::vec3(0.07f, beltY * 0.5f, beltHalfZ));

    double t = 0, spawnT = 0;
    auto step = [&]() {
        spawnT += 1.0 / 60;
        if (spawnT > 0.55 && (int)crates.size() < 40) { spawn(); spawnT = 0; }
        for (int s = 0; s < 4; s++) {
            wind.time = t;
            for (auto& o : crates) { belt.apply(o.body.get(), 1.0 / 240); wind.apply(o.body.get(), 1.0 / 240); }
            world.startFrame(); world.runPhysics(1.0 / 240); t += 1.0 / 240;
        }
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 1.0, 0), 5.5f);
        r.beginShadow();
        r.shadowDraw(box, beltM); r.shadowDraw(box, legL); r.shadowDraw(box, legR);
        for (auto& o : crates) r.shadowDraw(box, glm::scale(glmFromPhys(o.body->getTransform()), glm::vec3(o.half.x, o.half.y, o.half.z)));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.24, 0.23, 0.22), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, beltM, glm::vec3(0.14, 0.15, 0.17), 0.5f, 0.3f);
        r.drawPBR(r.pSolid, box, legL, glm::vec3(0.10, 0.10, 0.12), 0.6f, 0.4f);
        r.drawPBR(r.pSolid, box, legR, glm::vec3(0.10, 0.10, 0.12), 0.6f, 0.4f);
        for (auto& o : crates) r.drawPBR(r.pSolid, box, glm::scale(glmFromPhys(o.body->getTransform()), glm::vec3(o.half.x, o.half.y, o.half.z)), o.col, 0.6f, 0.0f);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame(); char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p); step(); }
        std::printf("wrote %d frames (crates=%d)\n", frames, (int)crates.size()); return 0;
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
