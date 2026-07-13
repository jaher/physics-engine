// Feature playground: everything from the framework-parity work in one scene —
// heightfield terrain, a raycast vehicle driving across it, a kinematic
// character controller walking and jumping, a motor-driven hinge windmill,
// capsules tumbling down the slopes, and a trigger zone that lights up when the
// car drives through. ./playground3d [--shot out.png [frames]] [--video f [n]]
#include "phys/phys.h"
#include "common/gfx.h"
#include <vector>
#include <cstring>
#include <cmath>
using namespace phys;

static glm::mat4 glmFromPhys(const Matrix4& t) {
    return glm::mat4(t.data[0], t.data[4], t.data[8], 0, t.data[1], t.data[5], t.data[9], 0,
                     t.data[2], t.data[6], t.data[10], 0, t.data[3], t.data[7], t.data[11], 1);
}

int main(int argc, char** argv) {
    bool headless = false; const char* shot = nullptr; const char* video = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--shot")) { headless = true; shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--video")) { headless = true; video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;

    // ---- heightfield terrain: gentle rolling dunes
    HeightField hf; hf.init(48, 48, 1.0, Vector3(-24, 0, -24));
    for (int z = 0; z < 48; z++) for (int x = 0; x < 48; x++)
        hf.at(x, z) = 0.5 * std::sin(x * 0.30) * std::cos(z * 0.24) + 0.25 * std::sin((x + z) * 0.5);

    // ---- raycast vehicle
    RigidBody chassis; chassis.setMass(1200);
    { Matrix3 I; I.setBlockInertiaTensor(Vector3(1.0, 0.4, 2.2), 1200); chassis.setInertiaTensor(I); }
    chassis.setPosition(-9, 1.6, -9); chassis.setAcceleration(0, -9.81, 0);
    chassis.setDamping(0.995, 0.75); chassis.setCanSleep(false); chassis.calculateDerivedData();
    RaycastVehicle car; car.chassis = &chassis; car.terrain = &hf;
    car.addWheel(Vector3(-0.85, -0.1, 1.35), true, false);
    car.addWheel(Vector3(0.85, -0.1, 1.35), true, false);
    car.addWheel(Vector3(-0.85, -0.1, -1.35), false, true);
    car.addWheel(Vector3(0.85, -0.1, -1.35), false, true);
    car.engineForce = 2600;

    // ---- character controller walking the dunes
    CharacterController hero; hero.position = Vector3(6, 2.0, -10); hero.terrain = &hf;

    // ---- windmill: motorized hinge about a world anchor atop a post
    Vector3 postTop(9, 4.0, 6);
    RigidBody blades; blades.setMass(40);
    { Matrix3 I; I.setBlockInertiaTensor(Vector3(2.4, 2.4, 0.15), 40); blades.setInertiaTensor(I); }
    blades.setPosition(postTop); blades.setDamping(0.99, 0.99); blades.setCanSleep(false);
    blades.setAcceleration(0, 0, 0); blades.calculateDerivedData();
    HingeJoint mill; mill.set(&blades, Vector3(0, 0, 0), Vector3(0, 0, 1), nullptr, postTop, Vector3(0, 0, 1));
    mill.hasMotor = true; mill.motorTargetVel = 1.6;
    JointSolver joints; joints.joints.push_back(&mill);

    // ---- capsules tumbling on the terrain
    struct Cap { RigidBody body; CollisionCapsule shape; glm::vec3 colour; };
    std::vector<Cap> caps(5);
    for (int i = 0; i < 5; i++) { Cap& c = caps[i];
        c.body.setMass(6); Matrix3 I; I.setBlockInertiaTensor(Vector3(0.25, 0.55, 0.25), 6); c.body.setInertiaTensor(I);
        c.body.setPosition(-4.0 + i * 2.0, 3.5 + i * 0.4, 4.0 + (i % 2));
        c.body.setOrientation(Quaternion(0.9, 0.25 * i, 0.2, 0.1)); c.body.setAcceleration(0, -9.81, 0);
        c.body.setDamping(0.6, 0.5); c.body.setCanSleep(false); c.body.calculateDerivedData();
        c.shape.body = &c.body; c.shape.radius = 0.24; c.shape.halfHeight = 0.36;
        c.colour = glm::vec3(0.3f + 0.12f * i, 0.5f - 0.06f * i, 0.7f - 0.08f * i);
    }
    struct CapGen : ContactGenerator {
        std::vector<Cap>* caps; const HeightField* hf;
        unsigned addContact(Contact* c, unsigned limit) const override {
            CollisionData d; d.contactArray = c; d.reset(limit); d.friction = 0.5; d.restitution = 0.1;
            for (auto& cp : *caps) { if (!d.hasMoreContacts()) break;
                cp.shape.calculateInternals(); hf->capsuleContact(cp.shape, &d); }
            for (size_t i = 0; i < caps->size(); i++) for (size_t j = i + 1; j < caps->size(); j++) {
                if (!d.hasMoreContacts()) break;
                CapsuleDetector::capsuleAndCapsule((*caps)[i].shape, (*caps)[j].shape, &d); }
            return d.contactCount;
        }
    } capGen; capGen.caps = &caps; capGen.hf = &hf;
    World world(400);
    for (auto& c : caps) world.getRigidBodies().push_back(&c.body);
    world.getContactGenerators().push_back(&capGen);

    // ---- trigger zone the car drives through
    TriggerVolume gate; gate.centre = Vector3(0, 1.0, 0); gate.isSphere = true; gate.radius = 3.2;
    bool gateLit = false;
    gate.onEnter = [&](RigidBody*) { gateLit = true; };
    gate.onExit = [&](RigidBody*) { gateLit = false; };
    std::vector<RigidBody*> carBody{&chassis};

    gfx::App app(W, H, "feature playground", headless);
    gfx::Renderer r; r.init(W, H);
    gfx::Mesh box = gfx::makeBox(), sphere = gfx::makeSphere(), plane = gfx::makePlane(60, 1);
    // terrain mesh
    gfx::Mesh terrainMesh;
    { std::vector<float> v; std::vector<unsigned> idx;
      for (int z = 0; z < 48; z++) for (int x = 0; x < 48; x++) {
          real wx = -24 + x, wz = -24 + z; Vector3 n = hf.normal(wx, wz);
          v.insert(v.end(), {(float)wx, (float)hf.sample(wx, wz), (float)wz, (float)n.x, (float)n.y, (float)n.z}); }
      for (int z = 0; z < 47; z++) for (int x = 0; x < 47; x++) {
          unsigned a = z * 48 + x, b = a + 1, c = a + 48, d = c + 1;
          idx.insert(idx.end(), {a, c, b, b, c, d}); }
      terrainMesh.upload(v, idx); }

    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 1.5, 0); cam.dist = 20.0f; cam.yaw = 0.6f; cam.pitch = 0.42f;
    glm::mat4 proj = glm::perspective(glm::radians(45.0f), (float)W / H, 0.1f, 200.0f);
    glfwSetWindowUserPointer(app.win, &cam);
    glfwSetScrollCallback(app.win, [](GLFWwindow* w, double, double dy) { auto* c = (gfx::OrbitCamera*)glfwGetWindowUserPointer(w); c->dist *= (dy > 0 ? 0.9f : 1.1f); });

    double t = 0;
    auto stepPhysics = [&](double time) {
        car.steer = 0.28 * std::sin(time * 0.5);           // wander across the dunes
        car.update(1.0 / 120); chassis.integrate(1.0 / 120);
        Vector3 fwd(std::sin(time * 0.35), 0, std::cos(time * 0.35));
        hero.move(fwd * 2.2, 1.0 / 120);
        if (std::fmod(time, 3.0) < 0.01) hero.jump(4.5);   // periodic hop
        blades.integrate(1.0 / 120); joints.solve(1.0 / 120);
        world.startFrame(); world.runPhysics(1.0 / 120);
        gate.frame(carBody);
    };
    auto renderFrame = [&]() {
        r.setLightForScene(glm::vec3(0, 1.5, 0), 20.0f);
        glm::mat4 chassisM = glm::scale(glmFromPhys(chassis.getTransform()), glm::vec3(0.95f, 0.45f, 2.1f));
        glm::mat4 bladeM1 = glm::scale(glmFromPhys(blades.getTransform()), glm::vec3(2.4f, 0.18f, 0.1f));
        glm::mat4 bladeM2 = glm::scale(glmFromPhys(blades.getTransform()), glm::vec3(0.18f, 2.4f, 0.1f));
        glm::mat4 postM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(postTop.x, (postTop.y - hf.sample(postTop.x, postTop.z)) * 0.5f + hf.sample(postTop.x, postTop.z), postTop.z)), glm::vec3(0.16f, (postTop.y - hf.sample(postTop.x, postTop.z)) * 0.5f, 0.16f));
        r.beginShadow();
        r.shadowDraw(terrainMesh, glm::mat4(1));
        r.shadowDraw(box, chassisM); r.shadowDraw(box, bladeM1); r.shadowDraw(box, bladeM2); r.shadowDraw(box, postM);
        for (auto& c : caps) r.shadowDraw(sphere, glm::scale(glmFromPhys(c.body.getTransform()), glm::vec3(0.24f, 0.6f, 0.24f)));
        r.endShadow();
        r.beginScene(cam.view(), proj, cam.eye());
        r.drawPBR(r.pSolid, terrainMesh, glm::mat4(1), glm::vec3(0.45, 0.52, 0.33), 0.95f, 0.0f);
        r.drawPBR(r.pSolid, box, chassisM, glm::vec3(0.75, 0.2, 0.15), 0.4f, 0.3f);
        for (auto& w : car.wheels) { Vector3 wc = car.wheelCentre(w);
            r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), glm::vec3(wc.x, wc.y, wc.z)), glm::vec3((float)w.radius)), glm::vec3(0.08, 0.08, 0.09), 0.8f, 0.0f); }
        // character capsule
        glm::mat4 heroM = glm::scale(glm::translate(glm::mat4(1), glm::vec3(hero.position.x, hero.position.y, hero.position.z)),
                                     glm::vec3((float)hero.radius, (float)(hero.halfHeight + hero.radius), (float)hero.radius));
        r.drawPBR(r.pSolid, sphere, heroM, glm::vec3(0.9, 0.6, 0.15), 0.5f, 0.0f);
        r.drawPBR(r.pSolid, box, postM, glm::vec3(0.4, 0.3, 0.2), 0.8f, 0.0f);
        r.drawPBR(r.pSolid, box, bladeM1, glm::vec3(0.85, 0.85, 0.9), 0.5f, 0.2f);
        r.drawPBR(r.pSolid, box, bladeM2, glm::vec3(0.85, 0.85, 0.9), 0.5f, 0.2f);
        for (auto& c : caps) r.drawPBR(r.pSolid, sphere,
            glm::scale(glmFromPhys(c.body.getTransform()), glm::vec3(0.24f, 0.6f, 0.24f)), c.colour, 0.5f, 0.0f);
        // trigger gate ring (lit when the car is inside)
        glm::vec3 gateCol = gateLit ? glm::vec3(0.2, 0.95, 0.3) : glm::vec3(0.5, 0.5, 0.55);
        for (int i = 0; i < 14; i++) { double a = i * 2 * M_PI / 14;
            glm::vec3 p(gate.centre.x + std::cos(a) * gate.radius, hf.sample(gate.centre.x + std::cos(a) * gate.radius, gate.centre.z + std::sin(a) * gate.radius) + 0.35, gate.centre.z + std::sin(a) * gate.radius);
            r.drawPBR(r.pSolid, sphere, glm::scale(glm::translate(glm::mat4(1), p), glm::vec3(0.16f)), gateCol, 0.4f, gateLit ? 0.4f : 0.0f); }
        r.endScene();
    };

    if (video) {
        for (int f = 0; f < frames; f++) { cam.yaw = 0.5f + f * (0.6f / frames);
            renderFrame(); char p[512]; std::snprintf(p, sizeof(p), "%s_%04d.png", video, f); r.screenshot(p);
            for (int s = 0; s < 2; s++) { stepPhysics(t); t += 1.0 / 120; } }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    if (headless) { for (int i = 0; i < frames * 2; i++) { stepPhysics(t); t += 1.0 / 120; } renderFrame(); r.screenshot(shot); std::printf("wrote %s\n", shot); return 0; }
    double lx = 0, ly = 0; bool drag = false;
    while (app.running()) {
        if (glfwGetKey(app.win, GLFW_KEY_ESCAPE) == GLFW_PRESS) break;
        double mx, my; glfwGetCursorPos(app.win, &mx, &my);
        bool down = glfwGetMouseButton(app.win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        if (down && drag) { cam.yaw -= (float)(mx - lx) * 0.005f; cam.pitch = glm::clamp(cam.pitch + (float)(my - ly) * 0.005f, -1.4f, 1.4f); }
        lx = mx; ly = my; drag = down;
        stepPhysics(t); t += 1.0 / 120; stepPhysics(t); t += 1.0 / 120;
        renderFrame(); r.present(); app.poll();
    }
    return 0;
}
