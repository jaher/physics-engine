// Galton board (bean machine / quincunx). Balls stream down through a triangular
// lattice of fixed pegs; at each peg a ball deflects left or right, and the many
// independent 50/50 choices make the balls pile into the bins below as a binomial
// → the bell curve of the central limit theorem. Everything runs on the engine's
// rigid-body dynamics and collision resolver (phys::World): balls are spheres,
// pegs are immovable spheres, the bins are immovable boxes, and the case is a set
// of half-space planes.
//   ./galton3d [--shot out.png [frames]] [--video f [n]]
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
static glm::vec3 hsv(float h, float s, float v) {
    float c = v * s, x = c * (1 - std::fabs(std::fmod(h * 6, 2.0f) - 1)), m = v - c;
    float r, g, b; int i = (int)(h * 6) % 6;
    if (i == 0) { r = c; g = x; b = 0; } else if (i == 1) { r = x; g = c; b = 0; }
    else if (i == 2) { r = 0; g = c; b = x; } else if (i == 3) { r = 0; g = x; b = c; }
    else if (i == 4) { r = x; g = 0; b = c; } else { r = c; g = 0; b = x; }
    return glm::vec3(r + m, g + m, b + m);
}

struct Fixed { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionSphere> csph; std::unique_ptr<CollisionBox> cbox; Vector3 half; real radius; };
struct Ball { std::unique_ptr<RigidBody> body; std::unique_ptr<CollisionSphere> csph; glm::vec3 col; int lastRow = 9999; };

struct Board : ContactGenerator {
    std::vector<Ball>* balls; std::vector<Fixed>* pegs; std::vector<Fixed>* divs;
    CollisionPlane floor, left, right, front, back; real fr, re;
    unsigned addContact(Contact* c, unsigned limit) const override {
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
        for (auto& b : *balls) b.csph->calculateInternals();
        auto& B = *balls;
        for (auto& b : B) {
            if (!d.hasMoreContacts()) break;
            CollisionDetector::sphereAndHalfSpace(*b.csph, floor, &d);
            CollisionDetector::sphereAndHalfSpace(*b.csph, left, &d);
            CollisionDetector::sphereAndHalfSpace(*b.csph, right, &d);
            CollisionDetector::sphereAndHalfSpace(*b.csph, front, &d);
            CollisionDetector::sphereAndHalfSpace(*b.csph, back, &d);
            for (auto& p : *pegs) { if (!d.hasMoreContacts()) break; CollisionDetector::sphereAndSphere(*b.csph, *p.csph, &d); }
            for (auto& v : *divs) { if (!d.hasMoreContacts()) break; CollisionDetector::boxAndSphere(*v.cbox, *b.csph, &d); }
        }
        for (size_t i = 0; i < B.size(); i++) for (size_t j = i + 1; j < B.size(); j++) {
            if (!d.hasMoreContacts()) break;
            if (std::fabs(B[i].body->getPosition().y - B[j].body->getPosition().y) > 0.3) continue;   // cheap y-band cull
            CollisionDetector::sphereAndSphere(*B[i].csph, *B[j].csph, &d);
        }
        return d.contactCount;
    }
};

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 700;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    const int ROWS = 11, NBALLS = 170;
    const real dx = 0.44, dy = 0.36, pr = 0.05, br = 0.075, zh = 0.12;
    const real binTop = 2.7, pegY0 = 3.05;                 // bottom peg row height
    const real half = 6.0 * dx + 0.11;                     // case half-width

    World world(7000);
    std::vector<Fixed> pegs, divs; std::vector<Ball> balls;

    auto fixSphere = [&](Vector3 pos, real rad) { Fixed f; f.radius = rad;
        f.body = std::make_unique<RigidBody>(); f.body->setInverseMass(0);
        Matrix3 I; I.setInertiaTensorCoeffs(1e12, 1e12, 1e12); f.body->setInertiaTensor(I);
        f.body->setPosition(pos); f.body->setDamping(1, 1); f.body->calculateDerivedData();
        f.csph = std::make_unique<CollisionSphere>(); f.csph->body = f.body.get(); f.csph->radius = rad; f.csph->calculateInternals();
        pegs.push_back(std::move(f)); };
    auto fixBox = [&](Vector3 pos, Vector3 hs) { Fixed f; f.half = hs;
        f.body = std::make_unique<RigidBody>(); f.body->setInverseMass(0);
        Matrix3 I; I.setInertiaTensorCoeffs(1e12, 1e12, 1e12); f.body->setInertiaTensor(I);
        f.body->setPosition(pos); f.body->setDamping(1, 1); f.body->calculateDerivedData();
        f.cbox = std::make_unique<CollisionBox>(); f.cbox->body = f.body.get(); f.cbox->halfSize = hs; f.cbox->calculateInternals();
        divs.push_back(std::move(f)); };

    // triangular quincunx of pegs (row r has r+1 pegs), rows counted from the top
    for (int r = 0; r < ROWS; r++) {
        real y = pegY0 + (ROWS - 1 - r) * dy;
        for (int j = 0; j <= r; j++) fixSphere(Vector3((j - r * 0.5) * dx, y, 0), pr);
    }
    // 13 bin dividers under the bottom row → 12 bins
    for (int i = -6; i <= 6; i++) fixBox(Vector3(i * dx, binTop * 0.5, 0), Vector3(0.014, binTop * 0.5, zh));

    Board board; board.balls = &balls; board.pegs = &pegs; board.divs = &divs;
    board.floor.direction = Vector3(0, 1, 0); board.floor.offset = 0;
    board.left.direction = Vector3(1, 0, 0);  board.left.offset = -half;
    board.right.direction = Vector3(-1, 0, 0); board.right.offset = -half;
    board.back.direction = Vector3(0, 0, 1);  board.back.offset = -zh;
    board.front.direction = Vector3(0, 0, -1); board.front.offset = -zh;
    board.fr = 0.5; board.re = 0.1;
    world.getContactGenerators().push_back(&board);

    unsigned sd = 12345u; auto rnd = [&]() { sd = sd * 1103515245u + 12345u; return (float)(((sd >> 16) & 0x7fff) / 32767.0f); };
    int spawned = 0; double spawnT = 0;
    auto spawn = [&]() {
        Ball b; b.col = hsv((float)spawned / NBALLS, 0.65f, 0.95f);
        b.body = std::make_unique<RigidBody>(); real m = br * br * br * 30; b.body->setMass(m);
        Matrix3 I; real c = 0.4 * m * br * br; I.setInertiaTensorCoeffs(c, c, c); b.body->setInertiaTensor(I);
        b.body->setPosition(Vector3((rnd() - 0.5f) * 0.02, pegY0 + (ROWS - 1) * dy + 0.55, (rnd() - 0.5f) * 0.04));
        b.body->setVelocity(0, 0, 0);
        b.body->setAcceleration(0, -9.81, 0); b.body->setDamping(0.995, 0.9); b.body->calculateDerivedData();
        b.csph = std::make_unique<CollisionSphere>(); b.csph->radius = br;
        balls.push_back(std::move(b)); balls.back().csph->body = balls.back().body.get();
        world.getRigidBodies().push_back(balls.back().body.get());
        spawned++;
    };

    gfx::App app(W, H, "galton board", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh sphere = gfx::makeSphere(), box = gfx::makeBox(), plane = gfx::makePlane(40, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 3.5, 0); cam.dist = 9.9f; cam.yaw = 0.09f; cam.pitch = 0.01f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy2) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy2 > 0 ? 0.9f : 1.1f); });

    // back board + side frame (visual case)
    glm::mat4 backM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(0, 3.4f, -zh - 0.03f)), glm::vec3(half + 0.12f, 3.7f, 0.03f));
    glm::mat4 frameL = glm::scale(glm::translate(glm::mat4(1), glm::vec3(-half - 0.05f, 3.4f, 0)), glm::vec3(0.05f, 3.7f, zh + 0.02f));
    glm::mat4 frameR = glm::scale(glm::translate(glm::mat4(1), glm::vec3(half + 0.05f, 3.4f, 0)), glm::vec3(0.05f, 3.7f, zh + 0.02f));

    const real topY = pegY0 + (ROWS - 1) * dy;
    auto stepPhysics = [&]() {
        spawnT += 1.0 / 60;
        if (spawned < NBALLS && spawnT > 0.05) { spawn(); spawnT = 0; }
        for (int s = 0; s < 4; s++) { world.startFrame(); world.runPhysics(1.0 / 240); }
        // micro-deflection as each ball reaches a new peg row → decorrelated 50/50
        // hops (models the imperfections a real board relies on) → a clean binomial
        for (auto& b : balls) {
            real y = b.body->getPosition().y;
            if (y <= binTop || y >= topY + 0.25) continue;
            int row = (int)std::floor((topY + 0.2 - y) / dy);
            if (row != b.lastRow) {                                  // fresh 50/50 hop of ~half a peg spacing
                Vector3 v = b.body->getVelocity();
                real hop = 0.5 * dx * std::max((real)1.4, std::fabs(v.y)) / dy;
                v.x = (rnd() < 0.5f ? -1 : 1) * hop; b.body->setVelocity(v); b.lastRow = row;
            }
        }
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 3.2, 0), 6.5f);
        r.beginShadow();
        r.shadowDraw(box, backM);
        for (auto& v : *(&divs)) { glm::mat4 m = glm::scale(glmFromPhys(v.body->getTransform()), glm::vec3(v.half.x, v.half.y, v.half.z)); r.shadowDraw(box, m); }
        for (auto& b : balls) { glm::mat4 m = glm::scale(glmFromPhys(b.body->getTransform()), glm::vec3(br)); r.shadowDraw(sphere, m); }
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.22, 0.20, 0.18), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, backM, glm::vec3(0.32, 0.26, 0.20), 0.85f, 0.0f);
        r.drawPBR(r.pSolid, box, frameL, glm::vec3(0.28, 0.22, 0.16), 0.8f, 0.0f);
        r.drawPBR(r.pSolid, box, frameR, glm::vec3(0.28, 0.22, 0.16), 0.8f, 0.0f);
        for (auto& v : divs) { glm::mat4 m = glm::scale(glmFromPhys(v.body->getTransform()), glm::vec3(v.half.x, v.half.y, v.half.z)); r.drawPBR(r.pSolid, box, m, glm::vec3(0.20, 0.17, 0.13), 0.7f, 0.0f); }
        for (auto& p : pegs) { glm::mat4 m = glm::scale(glmFromPhys(p.body->getTransform()), glm::vec3(pr)); r.drawPBR(r.pSolid, sphere, m, glm::vec3(0.55, 0.57, 0.6), 0.3f, 0.9f); }
        for (auto& b : balls) { glm::mat4 m = glm::scale(glmFromPhys(b.body->getTransform()), glm::vec3(br)); r.drawPBR(r.pSolid, sphere, m, b.col, 0.25f, 0.1f); }
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame(); char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p); stepPhysics(); }
        std::printf("wrote %d frames (balls=%d)\n", frames, (int)balls.size()); return 0;
    }
    if (headless) { for (int i = 0; i < frames; i++) stepPhysics(); renderFrame(); r.screenshot(shot);
        int hist[13] = {0}; for (auto& b : balls) { int bin = (int)std::floor(b.body->getPosition().x / dx + 6.5); if (bin >= 0 && bin < 13) hist[bin]++; }
        std::printf("wrote %s (balls=%d) bins:", shot, (int)balls.size()); for (int k = 0; k < 12; k++) std::printf(" %d", hist[k]); std::printf("\n"); return 0; }
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
