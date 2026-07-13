// A metal ball smashes THROUGH a concrete wall, leaving a jagged hole. The wall
// is pre-fractured into welded fragments (phys::fracture); on impact only the
// fragments inside a jagged, angle-modulated radius around the impact point are
// released and blasted forward — the rest of the wall stays intact. The ball
// punches through (losing some momentum) and the debris settles via the engine.
//   ./wall3d                 interactive        ./wall3d --shot out.png [frames]
//   ./wall3d --video f [n]
#include "phys/phys.h"
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
struct Frag { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionBox> box; Vector3 half; glm::vec3 colour; bool released = false; };

struct Collider : ContactGenerator {
    std::vector<Frag>* frags; RigidBody* ball; CollisionSphere* ballC; CollisionPlane ground;
    const bool* frozen; const bool* punched;
    real fr, re;
    unsigned addContact(Contact* c, unsigned limit) const override {
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
        ballC->calculateInternals();
        if (d.hasMoreContacts()) CollisionDetector::sphereAndHalfSpace(*ballC, ground, &d);
        if (*punched && !*frozen) {
            std::vector<Frag*> rel;
            for (auto& f : *frags) if (f.released) { f.box->calculateInternals(); rel.push_back(&f); }
            for (auto* f : rel) { if (!d.hasMoreContacts()) break; CollisionDetector::boxAndHalfSpace(*f->box, ground, &d); }
            for (size_t i = 0; i < rel.size(); i++) { if (!d.hasMoreContacts()) break;
                CollisionDetector::boxAndSphere(*rel[i]->box, *ballC, &d);
                for (size_t j = i + 1; j < rel.size(); j++) { if (!d.hasMoreContacts()) break;
                    CollisionDetector::boxAndBox(*rel[i]->box, *rel[j]->box, &d); } }
        }
        return d.contactCount;
    }
};

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 240;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    RigidBody::sleepEpsilon = (real)0.8;

    // ---- the wall: 4.4 m wide, 3 m tall, 0.24 m thick, fractured into a grid
    Vector3 wallC(0, 1.5, 0), wallHalf(0.85, 1.25, 0.018);   // a tall thin mirror
    auto descs = fractureBoxGrid(wallC, wallHalf, 24, 16, 1, 0.5, 31u);
    std::vector<Frag> frags;
    unsigned cs = 777; auto cr = [&]() { cs = cs * 1103515245u + 12345u; return (float)(((cs >> 16) & 0x7fff) / 32767.0); };
    for (auto& d : descs) {
        Frag f; f.half = d.halfSize;
        real vol = 8 * d.halfSize.x * d.halfSize.y * d.halfSize.z;
        f.body = std::make_unique<RigidBody>(); f.body->setMass(vol * 2.5);        // glass density
        Matrix3 I; I.setBlockInertiaTensor(d.halfSize, f.body->getMass()); f.body->setInertiaTensor(I);
        f.body->setPosition(d.centre); f.body->setAcceleration(0, -9.81, 0); f.body->setDamping(0.35, 0.3);
        f.body->calculateDerivedData(); f.body->setAwake(false);                    // welded into the wall
        f.box = std::make_unique<CollisionBox>(); f.box->body = f.body.get(); f.box->halfSize = d.halfSize;
        float j = (cr() - 0.5f) * 0.06f;
        f.colour = glm::vec3(0.78f + j, 0.82f + j, 0.88f + j);   // silvered glass
        frags.push_back(std::move(f));
    }

    // ---- the metal ball
    RigidBody ball; ball.setMass(8); real rr = 0.14;
    Matrix3 Ib; real cI = 0.4 * 8 * rr * rr; Ib.setInertiaTensorCoeffs(cI, cI, cI); ball.setInertiaTensor(Ib);
    ball.setAcceleration(0, -1.2, 0); ball.setDamping(0.999, 0.995);
    ball.setPosition(0.15, 1.7, -6.0); ball.setVelocity(0, 0.4, 14.0);             // flat shot at the wall centre
    ball.calculateDerivedData();
    CollisionSphere ballC; ballC.body = &ball; ballC.radius = rr;

    World world(5000);
    for (auto& f : frags) world.getRigidBodies().push_back(f.body.get());
    world.getRigidBodies().push_back(&ball);
    bool frozen = false, punched = false;
    Collider col; col.frags = &frags; col.ball = &ball; col.ballC = &ballC; col.frozen = &frozen; col.punched = &punched;
    col.ground.direction = Vector3(0, 1, 0); col.ground.offset = 0; col.fr = 0.75; col.re = 0.0;
    world.getContactGenerators().push_back(&col);

    gfx::App app(W, H, "mirror shatter", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.3, 0); cam.dist = 6.0f; cam.yaw = 2.85f; cam.pitch = 0.14f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    auto stepPhysics = [&]() {
        world.startFrame(); world.runPhysics(1.0 / 120);
        if (!punched && ball.getPosition().z > -wallHalf.z - rr - 0.02) {
            punched = true;
            Vector3 impact = ball.getPosition(); impact.z = wallC.z;
            // jagged hole: angle-modulated radius (three lobes + jitter per fragment)
            unsigned hs = 555; auto hr = [&]() { hs = hs * 1103515245u + 12345u; return (real)(((hs >> 16) & 0x7fff) / 32767.0); };
            real R = 9.0; real ph1 = hr() * 6.28, ph2 = hr() * 6.28;  // whole mirror lets go real ph1 = hr() * 6.28, ph2 = hr() * 6.28;
            Destructible dest;
            for (auto& f : frags) {
                Vector3 d = f.body->getPosition() - impact; d.z = 0;
                real ang = std::atan2(d.y, d.x);
                real rag = R * (1 + (real)0.26 * real_sin(3 * ang + ph1) + (real)0.14 * real_sin(5 * ang + ph2))
                         + (hr() - (real)0.5) * (real)0.14;
                if (d.magnitude() < rag) { f.released = true; dest.fragments.push_back(f.body.get()); }
            }
            dest.shatter(impact - Vector3(0, 0, 0.3), Vector3(0, 0.2, 1.2), 1.1, 10, 21u);
            ball.setVelocity(ball.getVelocity() * (real)0.55);        // punching through costs momentum
        } else if (punched && !frozen) {
            real ke = 0; int n = 0;
            for (auto& f : frags) if (f.released) { ke += f.body->getKineticEnergy(); n++; }
            if (n && ke < 0.05) { frozen = true; for (auto& f : frags) if (f.released) f.body->setAwake(false); }
        }
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 1.4, 0), 8.5f);
        r.beginShadow();
        for (auto& f : frags) { glm::mat4 m = glm::scale(glmFromPhys(f.body->getTransform()), glm::vec3(f.half.x, f.half.y, f.half.z)); r.shadowDraw(box, m); }
        r.shadowDraw(sphere, glm::scale(glmFromPhys(ball.getTransform()), glm::vec3(rr)));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.56, 0.56, 0.53), 0.95f, 0.0f);
        for (auto& f : frags) { glm::mat4 m = glm::scale(glmFromPhys(f.body->getTransform()), glm::vec3(f.half.x, f.half.y, f.half.z));
            r.drawPBR(r.pSolid, box, m, f.colour, 0.06f, 0.95f); }   // mirror finish
        r.drawPBR(r.pSolid, sphere, glm::scale(glmFromPhys(ball.getTransform()), glm::vec3(rr)), glm::vec3(0.62, 0.64, 0.68), 0.25f, 1.0f);
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 2.75f - f * (0.5f / frames); renderFrame();
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
