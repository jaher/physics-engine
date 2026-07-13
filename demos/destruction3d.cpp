// Object destruction: a solid block (granite or wood) shatters when a heavy ball
// smashes into it. The block is pre-split into welded box fragments (phys::
// fracture); on impact they burst apart and then behave as ordinary rigid bodies,
// colliding and settling via the engine. Granite → chunky shards; wood → long
// splinters along the grain.
//   ./destruction3d [--material granite|wood]
//   ./destruction3d --material wood --shot out.png [frames]   |   --video f [n]
#include "phys/phys.h"
#include "common/gfx.h"
#include <vector>
#include <memory>
#include <cstring>
using namespace phys;

static glm::mat4 glmFromPhys(const Matrix4& t) {
    return glm::mat4(t.data[0], t.data[4], t.data[8], 0, t.data[1], t.data[5], t.data[9], 0,
                     t.data[2], t.data[6], t.data[10], 0, t.data[3], t.data[7], t.data[11], 1);
}
struct Frag { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionBox> box; Vector3 half; glm::vec3 colour; };

struct Collider : ContactGenerator {
    std::vector<Frag>* frags; RigidBody* ball; CollisionSphere* ballC; CollisionPlane ground;
    Destructible* dest; real fr, re; const bool* frozen = nullptr;
    unsigned addContact(Contact* c, unsigned limit) const override {
        if (!dest->shattered || (frozen && *frozen)) return 0;
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
        for (auto& f : *frags) f.box->calculateInternals();
        ballC->calculateInternals();
        for (auto& f : *frags) { if (!d.hasMoreContacts()) break; CollisionDetector::boxAndHalfSpace(*f.box, ground, &d); }
        if (d.hasMoreContacts()) CollisionDetector::sphereAndHalfSpace(*ballC, ground, &d);
        for (size_t i = 0; i < frags->size(); i++) { if (!d.hasMoreContacts()) break;
            CollisionDetector::boxAndSphere(*(*frags)[i].box, *ballC, &d);
            for (size_t j = i + 1; j < frags->size(); j++) { if (!d.hasMoreContacts()) break;
                CollisionDetector::boxAndBox(*(*frags)[i].box, *(*frags)[j].box, &d); }
        }
        return d.contactCount;
    }
};

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 260;
    std::string material = "granite";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--material")) material = argv[++i];
        else if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    bool wood = (material == "wood");
    const int W = 1280, H = 800;
    // Let settled debris deactivate readily so it doesn't micro-jitter forever
    // (box-box gives one contact per pair, so resting stacks jiggle otherwise).
    // Long wood splinters rock more, so they sleep at a slightly higher threshold.
    RigidBody::sleepEpsilon = wood ? (real)0.9 : (real)0.7;

    Vector3 blockC(0, 0, 0), blockHalf;
    int nx, ny, nz; real jitter, density; glm::vec3 baseCol; real burst, spin;
    Vector3 approach;                                        // ball travel direction
    // wood: grain along x → few cuts along x, many across → long splinters
    if (wood) { blockHalf = Vector3(1.3, 0.36, 0.36); nx = 3; ny = 6; nz = 6; jitter = 0.55; density = 0.7;
                baseCol = glm::vec3(0.46, 0.29, 0.14); burst = 2.6; spin = 8; approach = Vector3(0, 0, 1); }
    else      { blockHalf = Vector3(0.72, 0.78, 0.72); nx = 4; ny = 5; nz = 4; jitter = 0.42; density = 2.6;
                baseCol = glm::vec3(0.52, 0.51, 0.53); burst = 3.0; spin = 8; approach = Vector3(1, 0, 0); }
    blockC.y = blockHalf.y + 0.02;                           // block resting on the ground

    auto descs = fractureBoxGrid(blockC, blockHalf, nx, ny, nz, jitter, wood ? 99u : 7u);
    std::vector<Frag> frags;
    Destructible dest;
    unsigned cs = 4242; auto cr = [&]() { cs = cs * 1103515245u + 12345u; return (float)(((cs >> 16) & 0x7fff) / 32767.0); };
    for (auto& d : descs) {
        Frag f; f.half = d.halfSize;
        real vol = 8 * d.halfSize.x * d.halfSize.y * d.halfSize.z;
        f.body = std::make_unique<RigidBody>(); f.body->setMass(vol * density);
        Matrix3 I; I.setBlockInertiaTensor(d.halfSize, f.body->getMass()); f.body->setInertiaTensor(I);
        f.body->setPosition(d.centre); f.body->setAcceleration(0, -9.81, 0);
        f.body->setDamping(wood ? 0.3 : 0.35, wood ? 0.1 : 0.35);       // heavy angular damping tames splinter rocking
        f.body->calculateDerivedData(); f.body->setAwake(false);        // welded until shatter
        f.box = std::make_unique<CollisionBox>(); f.box->body = f.body.get(); f.box->halfSize = d.halfSize;
        float j = (cr() - 0.5f) * (wood ? 0.14f : 0.12f), grain = wood ? (cr() - 0.5f) * 0.08f : 0;
        f.colour = baseCol + glm::vec3(j + grain, j, j - grain * 0.5f);
        frags.push_back(std::move(f));
    }
    for (auto& f : frags) dest.fragments.push_back(f.body.get());

    // the wrecking ball
    RigidBody ball; ball.setMass(60); real rr = 0.45; Matrix3 Ib; real cI = 0.4 * 60 * rr * rr; Ib.setInertiaTensorCoeffs(cI, cI, cI);
    ball.setInertiaTensor(Ib); ball.setAcceleration(0, -2.0, 0); ball.setDamping(0.999, 0.999);
    Vector3 start = blockC - approach * 6.0; ball.setPosition(start);
    ball.setVelocity(approach * 13.0 + Vector3(0, 0.6, 0)); ball.calculateDerivedData();
    CollisionSphere ballC; ballC.body = &ball; ballC.radius = rr;

    World world(4000);                                      // iterations = 0 → auto-scaled to the contact count
    for (auto& f : frags) world.getRigidBodies().push_back(f.body.get());
    world.getRigidBodies().push_back(&ball);
    bool frozen = false;
    Collider col; col.frags = &frags; col.ball = &ball; col.ballC = &ballC; col.dest = &dest; col.frozen = &frozen;
    col.ground.direction = Vector3(0, 1, 0); col.ground.offset = 0; col.fr = 0.72; col.re = 0.0;
    world.getContactGenerators().push_back(&col);

    gfx::App app(W, H, wood ? "wood destruction" : "granite destruction", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 0.7, 0); cam.dist = 6.8f; cam.yaw = wood ? 1.15f : 0.75f; cam.pitch = 0.2f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    real halfAlong = blockHalf.x * real_abs(approach.x) + blockHalf.y * real_abs(approach.y) + blockHalf.z * real_abs(approach.z);
    auto stepPhysics = [&]() {
        world.startFrame(); world.runPhysics(1.0 / 120);
        if (!dest.shattered) {
            Vector3 rel = ball.getPosition() - blockC;
            real along = rel.x * approach.x + rel.y * approach.y + rel.z * approach.z;
            if (along > -halfAlong - rr) {
                Vector3 impact = blockC - approach * halfAlong; impact.y = ball.getPosition().y;
                dest.shatter(impact, approach * (burst * 0.55) + Vector3(0, burst * 0.5, 0), burst, spin, wood ? 5u : 3u);
            }
        } else if (!frozen) {
            // Once the debris island has come to rest, freeze it: interlocked box
            // fragments keep waking each other, so per-body sleep never converges.
            real ke = 0; for (auto& f : frags) ke += f.body->getKineticEnergy();
            if (ke < 0.06) { frozen = true; for (auto& f : frags) f.body->setAwake(false); }
        }
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 0.7, 0), 6.5f);
        r.beginShadow();
        for (auto& f : frags) { glm::mat4 m = glm::scale(glmFromPhys(f.body->getTransform()), glm::vec3(f.half.x, f.half.y, f.half.z)); r.shadowDraw(box, m); }
        r.shadowDraw(sphere, glm::scale(glmFromPhys(ball.getTransform()), glm::vec3(rr)));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.55, 0.55, 0.52), 0.95f, 0.0f);
        for (auto& f : frags) { glm::mat4 m = glm::scale(glmFromPhys(f.body->getTransform()), glm::vec3(f.half.x, f.half.y, f.half.z));
            r.drawPBR(r.pSolid, box, m, f.colour, wood ? 0.8f : 0.7f, 0.0f); }
        r.drawPBR(r.pSolid, sphere, glm::scale(glmFromPhys(ball.getTransform()), glm::vec3(rr)), glm::vec3(0.08, 0.08, 0.09), 0.35f, 0.6f);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw += 0.35f / frames; renderFrame();
            char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 2; s++) stepPhysics(); }
        std::printf("wrote %d frames\n", frames); return 0;
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
