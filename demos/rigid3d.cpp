// 3D rigid-body demo: boxes and spheres drop and pile up on a floor, driven by
// the engine's rigid-body dynamics, collision detection and contact resolver.
// Rendered with PBR + shadow mapping.  Interactive (drag to orbit, scroll to
// zoom), or headless:  ./rigid3d --shot out.png [frames]
#include "phys/phys.h"
#include "common/gfx.h"
#include <vector>
#include <memory>
#include <cstring>
using namespace phys;

static glm::mat4 glmFromPhys(const Matrix4& t) {
    return glm::mat4(t.data[0], t.data[4], t.data[8], 0,
                     t.data[1], t.data[5], t.data[9], 0,
                     t.data[2], t.data[6], t.data[10], 0,
                     t.data[3], t.data[7], t.data[11], 1);
}

struct Obj { std::unique_ptr<RigidBody> body; bool box; Vector3 half; real radius; glm::vec3 colour;
             std::unique_ptr<CollisionBox> cbox; std::unique_ptr<CollisionSphere> csph; };

struct SceneCollider : ContactGenerator {
    std::vector<Obj>* objs; CollisionPlane ground; real fr, re;
    unsigned addContact(Contact* c, unsigned limit) const override {
        CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
        for (auto& o : *objs) { if (o.box) o.cbox->calculateInternals(); else o.csph->calculateInternals(); }
        for (auto& o : *objs) { if (!d.hasMoreContacts()) break;
            if (o.box) CollisionDetector::boxAndHalfSpace(*o.cbox, ground, &d);
            else CollisionDetector::sphereAndHalfSpace(*o.csph, ground, &d); }
        for (size_t i = 0; i < objs->size(); i++) for (size_t j = i + 1; j < objs->size(); j++) {
            if (!d.hasMoreContacts()) break; Obj& a = (*objs)[i]; Obj& b = (*objs)[j];
            if (a.box && b.box) CollisionDetector::boxAndBox(*a.cbox, *b.cbox, &d);
            else if (a.box && !b.box) CollisionDetector::boxAndSphere(*a.cbox, *b.csph, &d);
            else if (!a.box && b.box) CollisionDetector::boxAndSphere(*b.cbox, *a.csph, &d);
            else CollisionDetector::sphereAndSphere(*a.csph, *b.csph, &d);
        }
        return d.contactCount;
    }
};

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; int frames = 260;
    for (int i = 1; i < argc; i++) { if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); } }
    const int W = 1280, H = 800;

    World world(200);
    std::vector<Obj> objs;
    auto addBox = [&](Vector3 pos, Vector3 half, glm::vec3 col, Quaternion q = Quaternion()) {
        Obj o; o.box = true; o.half = half; o.colour = col;
        o.body = std::make_unique<RigidBody>(); o.body->setMass(half.x * half.y * half.z * 8);
        Matrix3 I; I.setBlockInertiaTensor(half, o.body->getMass()); o.body->setInertiaTensor(I);
        o.body->setPosition(pos); o.body->setOrientation(q); o.body->setAcceleration(0, -9.81, 0);
        o.body->setDamping(0.995, 0.9); o.body->calculateDerivedData();
        o.cbox = std::make_unique<CollisionBox>(); o.cbox->body = o.body.get(); o.cbox->halfSize = half;
        objs.push_back(std::move(o)); };
    auto addSphere = [&](Vector3 pos, real rad, glm::vec3 col) {
        Obj o; o.box = false; o.radius = rad; o.colour = col;
        o.body = std::make_unique<RigidBody>(); real m = rad * rad * rad * 30; o.body->setMass(m);
        Matrix3 I; real c = 0.4 * m * rad * rad; I.setInertiaTensorCoeffs(c, c, c); o.body->setInertiaTensor(I);
        o.body->setPosition(pos); o.body->setAcceleration(0, -9.81, 0);
        o.body->setDamping(0.995, 0.9); o.body->calculateDerivedData();
        o.csph = std::make_unique<CollisionSphere>(); o.csph->body = o.body.get(); o.csph->radius = rad;
        objs.push_back(std::move(o)); };

    // a scattered drop that piles up
    addBox(Vector3(0, 0.6, 0), Vector3(1.2, 0.4, 1.2), glm::vec3(0.55, 0.15, 0.12));
    addBox(Vector3(-0.3, 1.6, 0.2), Vector3(0.6, 0.6, 0.6), glm::vec3(0.15, 0.32, 0.6), Quaternion(0.96, 0.2, 0.1, 0));
    addBox(Vector3(0.5, 2.8, -0.2), Vector3(0.5, 0.5, 0.5), glm::vec3(0.85, 0.7, 0.2), Quaternion(0.9, 0.1, 0.3, 0.1));
    addBox(Vector3(-1.2, 3.6, 0.6), Vector3(0.45, 0.9, 0.45), glm::vec3(0.2, 0.5, 0.25));
    addBox(Vector3(1.0, 4.8, 0.4), Vector3(0.7, 0.35, 0.7), glm::vec3(0.5, 0.3, 0.55), Quaternion(0.8, 0.4, 0.2, 0.3));
    addSphere(Vector3(0.1, 3.9, 0.1), 0.55, glm::vec3(0.75, 0.4, 0.15));
    addSphere(Vector3(-0.8, 5.6, -0.4), 0.5, glm::vec3(0.3, 0.55, 0.62));
    addSphere(Vector3(1.4, 6.5, 0.2), 0.45, glm::vec3(0.7, 0.72, 0.75));
    addBox(Vector3(0.0, 7.5, 0.0), Vector3(0.6, 0.6, 0.6), glm::vec3(0.62, 0.6, 0.15), Quaternion(0.85, 0.3, 0.3, 0.2));

    SceneCollider collider; collider.objs = &objs; collider.ground.direction = Vector3(0, 1, 0); collider.ground.offset = 0;
    collider.fr = 0.6; collider.re = 0.15;
    for (auto& o : objs) world.getRigidBodies().push_back(o.body.get());
    world.getContactGenerators().push_back(&collider);

    gfx::App app(W, H, "rigid bodies", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(40, 1);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.2, 0); cam.dist = 9.5f; cam.yaw = 0.7f; cam.pitch = 0.28f;
    glm::mat4 proj = glm::perspective(glm::radians(42.0f), (float)W / H, 0.1f, 100.0f);

    // mouse orbit
    double lx = 0, ly = 0; bool dragging = false;
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) {
        auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });
    glfwSetWindowUserPointer(app.win, &cam);

    auto renderFrame = [&]() {
        glm::vec3 camPos = cam.eye();
        r.setLightForScene(glm::vec3(0, 1.5, 0), 9.0f);
        r.beginShadow();
        for (auto& o : objs) { glm::mat4 m = glmFromPhys(o.body->getTransform());
            m = glm::scale(m, o.box ? glm::vec3(o.half.x, o.half.y, o.half.z) : glm::vec3(o.radius));
            r.shadowDraw(o.box ? box : sphere, m); }
        r.endShadow();
        r.beginScene(cam.view(), proj, camPos);
        r.drawPBR(r.pGround, plane, glm::mat4(1), glm::vec3(0.62, 0.60, 0.56), 0.95f, 0.0f);
        for (auto& o : objs) { glm::mat4 m = glmFromPhys(o.body->getTransform());
            m = glm::scale(m, o.box ? glm::vec3(o.half.x, o.half.y, o.half.z) : glm::vec3(o.radius));
            r.drawPBR(r.pSolid, o.box ? box : sphere, m, o.colour, o.box ? 0.55f : 0.3f, 0.0f); }
        r.endScene();
    };

    if (headless) {
        for (int i = 0; i < frames; i++) { world.startFrame(); world.runPhysics(1.0 / 120); world.startFrame(); world.runPhysics(1.0 / 120); }
        renderFrame();
        r.screenshot(shot); std::printf("wrote %s\n", shot); return 0;
    }
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && dragging) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch += (float)(my - ly) * 0.005f;
            cam.pitch = glm::clamp(cam.pitch, -1.4f, 1.4f); }
        lx = mx; ly = my; dragging = down;
        world.startFrame(); world.runPhysics(1.0 / 120);
        renderFrame(); r.present(); app.poll();
    }
    return 0;
}
