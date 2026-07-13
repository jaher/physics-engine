// Tests for the soft-constraint / actuation module: CFM softness → less force,
// a breakable joint that fires exactly on crossing, a conveyor spinning a box up
// to belt speed, wind pushing a body downwind, a radial field's falloff, a dense
// body sinking to a bounded terminal velocity, and an implicit spring that stays
// stable at a huge k·dt² where explicit Euler blows up.
#include "check.h"
#include "phys/body.h"
#include "phys/constraint2.h"
#include <cmath>
using namespace phys;

static void setBody(RigidBody& b, Vector3 pos, real mass = 1) {
    b.setMass(mass);
    Matrix3 I; I.setBlockInertiaTensor(Vector3(0.5, 0.5, 0.5), mass); b.setInertiaTensor(I);
    b.setPosition(pos); b.setDamping(1, 1); b.setCanSleep(false); b.calculateDerivedData();
}

int main() {
    // ---- 1. soft distance constraint: higher CFM softness → smaller force ----
    {
        real softs[2] = {0, 5}, forces[2];
        for (int s = 0; s < 2; s++) {
            RigidBody b; setBody(b, Vector3(0, -1.2, 0));           // hangs 1.2 below anchor, rest 1.0
            SoftDistanceConstraint c; c.set(&b, Vector3(0, 0, 0), nullptr, Vector3(0, 0, 0), 1.0);
            c.softness = softs[s];
            c.solve(0.01);
            forces[s] = c.lastForce;
        }
        CHECK(forces[0] > forces[1]);                              // softer link carries less force
        CHECK(forces[1] > 0);                                      // ...but still restores
    }

    // ---- 2. breakable joint: fires exactly when load crosses breakForce ----
    {
        int fires = 0;
        RigidBody b; setBody(b, Vector3(0, -1.0, 0));
        SoftDistanceConstraint c; c.set(&b, Vector3(0, 0, 0), nullptr, Vector3(0, 0, 0), 1.0);
        c.softness = 0; c.erp = 0.2;
        BreakableJoint bj; bj.constraint = &c; bj.breakForce = 250;  // rigid link, load = 2000·stretch
        bj.onBreak = [&fires]() { fires++; };
        real dt = 0.01;
        // ramp the stretch e from 0.05 upward; load = 2000·e crosses 250 at e=0.125
        for (int i = 5; i <= 20; i++) {
            real e = i * (real)0.01;
            b.setVelocity(0, 0, 0); b.setPosition(0, -(1.0 + e), 0); b.calculateDerivedData();
            bool wasBroken = bj.broken;
            bj.solve(dt);
            if (!wasBroken && e < (real)0.124) CHECK(!bj.broken);   // must not break before crossing
            if (!wasBroken && e < (real)0.124) CHECK(fires == 0);
        }
        CHECK(bj.broken);                                          // broke once past the threshold
        CHECK(fires == 1);                                         // fired exactly once
        for (int i = 0; i < 10; i++) bj.solve(dt);                 // extra steps must not re-fire
        CHECK(fires == 1);

        int fires2 = 0;                                            // strong link never breaks
        RigidBody b2; setBody(b2, Vector3(0, -1.3, 0));
        SoftDistanceConstraint c2; c2.set(&b2, Vector3(0, 0, 0), nullptr, Vector3(0, 0, 0), 1.0);
        BreakableJoint bj2; bj2.constraint = &c2; bj2.breakForce = (real)1e9;
        bj2.onBreak = [&fires2]() { fires2++; };
        for (int i = 0; i < 50; i++) { b2.setVelocity(0, 0, 0); b2.setPosition(0, -1.3, 0); b2.calculateDerivedData(); bj2.solve(dt); }
        CHECK(!bj2.broken); CHECK(fires2 == 0);
    }

    // ---- 3. conveyor drives a resting box up to belt speed ----
    {
        RigidBody box; setBody(box, Vector3(0, 0.5, 0), 3.0);
        box.setAcceleration(0, 0, 0); box.setVelocity(0, 0, 0);
        ConveyorSurface belt; belt.normal = Vector3(0, 1, 0); belt.offset = 0;
        belt.beltVelocity = Vector3(2, 0, 0); belt.maxHeight = 1.0; belt.mu = 1.0;
        for (int i = 0; i < 400; i++) { belt.apply(&box, 0.01); box.integrate(0.01); }
        CHECK_NEAR(box.getVelocity().x, 2.0, 0.02);               // reached belt speed
        CHECK_NEAR(box.getVelocity().y, 0.0, 1e-9);               // no vertical drive
        CHECK(box.getPosition().x > 0.5);                          // and moved downstream
    }

    // ---- 4a. wind pushes a body downwind ----
    {
        RigidBody b; setBody(b, Vector3(0, 0, 0), 1.0); b.setAcceleration(0, 0, 0);
        WindField wind; wind.wind = Vector3(5, 0, 0); wind.gustAmp = 0;
        wind.linearDrag = 0.8; wind.quadraticDrag = 0.2;
        for (int i = 0; i < 300; i++) { wind.apply(&b, 0.01); b.integrate(0.01); }
        CHECK(b.getVelocity().x > 0);                             // accelerated downwind
        CHECK(b.getPosition().x > 0);                             // and carried downwind
        CHECK(b.getVelocity().x < 5.0 + 1e-9);                    // approaches, never exceeds wind speed
        WindField gust; gust.wind = Vector3(0, 0, 0); gust.gustAxis = Vector3(1, 0, 0);
        gust.gustAmp = 3; gust.gustFreq = 2 * real_pi;
        CHECK_NEAR(gust.velocityAt(0.25).x, 3.0, 1e-9);          // gust peaks a quarter period in
    }

    // ---- 4b. radial field: attract / repel + falloff ----
    {
        RigidBody b; setBody(b, Vector3(2, 0, 0), 1.0); b.setAcceleration(0, 0, 0);
        RadialField att; att.centre = Vector3(0, 0, 0); att.strength = 4;
        for (int i = 0; i < 50; i++) { att.apply(&b, 0.01); b.integrate(0.01); }
        CHECK(b.getVelocity().x < 0);                             // pulled toward the centre
        RigidBody r; setBody(r, Vector3(2, 0, 0), 1.0); r.setAcceleration(0, 0, 0);
        RadialField rep; rep.centre = Vector3(0, 0, 0); rep.strength = -4;
        for (int i = 0; i < 50; i++) { rep.apply(&r, 0.01); r.integrate(0.01); }
        CHECK(r.getVelocity().x > 0);                             // pushed away
        RadialField inv; inv.strength = 1; inv.inverseSquare = true;
        RadialField lin; lin.strength = 1; lin.inverseSquare = false;
        CHECK_NEAR(inv.forceAt(Vector3(4, 0, 0)).magnitude(), 1.0 / 16.0, 1e-9);  // 1/r²
        CHECK_NEAR(lin.forceAt(Vector3(4, 0, 0)).magnitude(), 4.0, 1e-9);         // ∝ r
        CHECK(lin.forceAt(Vector3(4, 0, 0)).magnitude() > inv.forceAt(Vector3(4, 0, 0)).magnitude());
    }

    // ---- 5. dense body sinks to a bounded terminal velocity ----
    {
        RigidBody b; setBody(b, Vector3(0, 0, 0), 5.0);           // 5 kg, 1 L → density 5000
        b.setAcceleration(0, -9.81, 0);
        RigidHydrodynamics water; water.fluidDensity = 1000; water.waterHeight = 100;
        water.volume = 0.001; water.bodyHeight = 0.1; water.dragCoeff = 1;
        water.addedMass = 0.5; water.liftCoeff = 0; water.g = Vector3(0, -9.81, 0);
        real dt = 0.005;
        for (int i = 0; i < 4000; i++) { water.apply(&b, dt); b.integrate(dt); }
        real vy = b.getVelocity().y;
        CHECK(std::isfinite(vy));                                 // did not diverge
        CHECK(vy < 0);                                            // sinking
        CHECK(vy > -1.0);                                         // bounded terminal speed
        real vterm = -std::sqrt(((5.0 - 1.0) * 9.81) / (0.5 * 1000.0 * 1.0));
        CHECK_NEAR(vy, vterm, 0.02);                              // matches analytic terminal velocity
    }

    // ---- 6. implicit spring stays bounded where explicit Euler blows up ----
    {
        real dt = 0.1, k = 1e4, m = 1;                            // k·dt² = 100 ≫ 4 → explicit unstable
        ImplicitSpring1D imp; imp.mass = m; imp.k = k; imp.c = 0; imp.x = 1; imp.v = 0;
        ExplicitSpring1D exp_; exp_.mass = m; exp_.k = k; exp_.c = 0; exp_.x = 1; exp_.v = 0;
        real E0 = imp.energy();
        for (int i = 0; i < 200; i++) { imp.step(dt); exp_.step(dt); }
        CHECK(std::isfinite(imp.energy()));                       // implicit: finite
        CHECK(imp.energy() <= E0 + 1e-6);                         // ...bounded (numerically dissipative)
        CHECK(std::fabs(imp.x) < 2.0);                            // ...position stays put
        CHECK(exp_.energy() > imp.energy() * 1e3);               // explicit: blew up
        // small chain: same stability at large stiff dt
        ImplicitSpringChain chain; chain.mass = 1; chain.k = 1e4; chain.c = 0; chain.resize(5);
        chain.u = {0.5, -0.3, 0.2, -0.1, 0.4};
        real Ec0 = chain.energy();
        for (int i = 0; i < 300; i++) chain.step(0.1);
        CHECK(std::isfinite(chain.energy()));
        CHECK(chain.energy() <= Ec0 + 1e-6);                     // whole chain stays bounded
    }

    return test::report("constraint2");
}
