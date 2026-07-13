#include "phys/pworld.h"
#include "check.h"
using namespace phys;

int main() {
    // A) Projectile: semi-implicit Euler gives exact velocity and near-exact position.
    {
        Particle p; p.setMass(1); p.setDamping(1); p.setAcceleration(0, -9.81, 0); p.setVelocity(10, 20, 0);
        real dt = 0.001; for (int i = 0; i < 1000; i++) p.integrate(dt);   // t = 1s
        CHECK_NEAR(p.getVelocity().x, 10.0, 1e-9);
        CHECK_NEAR(p.getVelocity().y, 20.0 - 9.81, 1e-6);
        CHECK_NEAR(p.getPosition().x, 10.0, 1e-6);
        CHECK_NEAR(p.getPosition().y, 20.0 - 0.5 * 9.81, 2e-2);            // 15.095
    }

    // B) A mass hanging on an anchored spring settles at len = rest + mg/k.
    {
        Vector3 anchor(0, 0, 0);
        ParticleAnchoredSpring spring(&anchor, 20.0, 1.0);
        Particle p; p.setMass(2); p.setDamping(0.4); p.setAcceleration(0, -9.81, 0); p.setPosition(0, -1.0, 0);
        ParticleWorld world(1);
        world.getParticles().push_back(&p);
        world.getForceRegistry().add(&p, &spring);
        for (int i = 0; i < 8000; i++) { world.startFrame(); world.runPhysics(0.01); }
        CHECK_NEAR(-p.getPosition().y, 1.0 + 2.0 * 9.81 / 20.0, 2e-2);     // 1.981
    }

    // C) A rod keeps two particles a fixed distance apart.
    {
        Particle a, b;
        a.setMass(1); b.setMass(1); a.setDamping(1); b.setDamping(1);
        a.setPosition(0, 0, 0); b.setPosition(2, 0, 0);
        b.setVelocity(0, 5, 0);
        ParticleRod rod; rod.particle[0] = &a; rod.particle[1] = &b; rod.length = 2;
        ParticleWorld world(10);
        world.getParticles().push_back(&a); world.getParticles().push_back(&b);
        world.getContactGenerators().push_back(&rod);
        real maxErr = 0;
        for (int i = 0; i < 400; i++) { world.startFrame(); world.runPhysics(0.01);
            maxErr = std::max(maxErr, std::fabs((a.getPosition() - b.getPosition()).magnitude() - 2.0)); }
        CHECK(maxErr < 0.05);
    }

    // D) A ball dropped from h=5 with restitution 0.8 rebounds to ~e^2 h = 3.2.
    {
        Particle p; p.setMass(1); p.setDamping(1); p.setAcceleration(0, -9.81, 0); p.setPosition(0, 5, 0);
        ParticleWorld world(1);
        world.getParticles().push_back(&p);
        GroundContacts ground; ground.restitution = 0.8; ground.init(&world.getParticles());
        world.getContactGenerators().push_back(&ground);
        bool bounced = false; real peak = 0;
        for (int i = 0; i < 4000; i++) { world.startFrame(); world.runPhysics(0.002);
            if (!bounced && p.getPosition().y < 0.05 && p.getVelocity().y > 0) bounced = true;
            if (bounced) peak = std::max(peak, p.getPosition().y); }
        CHECK_NEAR(peak, 0.64 * 5.0, 0.4);
    }

    // E) A cable never lets the particle exceed its length from the anchor.
    {
        Particle p; p.setMass(1); p.setDamping(1); p.setAcceleration(0, -9.81, 0);
        p.setPosition(0, -2, 0); p.setVelocity(6, 0, 0);
        ParticleCableConstraint cable; cable.particle = &p; cable.anchor = Vector3(0, 0, 0);
        cable.maxLength = 2; cable.restitution = 0;
        ParticleWorld world(1);
        world.getParticles().push_back(&p);
        world.getContactGenerators().push_back(&cable);
        real maxDist = 0;
        for (int i = 0; i < 2000; i++) { world.startFrame(); world.runPhysics(0.005);
            maxDist = std::max(maxDist, p.getPosition().magnitude()); }
        CHECK(maxDist < 2.05);
    }

    return test::report("particle");
}
