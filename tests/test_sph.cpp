// Non-Newtonian SPH: the shear-rate-dependent rheology has the right shape, a
// pool is stable/conservative, and a shear-thickening fluid resists an impactor
// more than a Newtonian one.
#include "phys/sph.h"
#include "check.h"
#include <algorithm>
#include <cmath>
using namespace phys;

// Drop an identical ball into a small pool; return where it comes to rest.
static double dropRestY(bool water) {
    SPHFluid g; g.rheo = water ? Rheology::water() : Rheology::oobleck();
    Vector3 cmin(-0.3, 0, -0.2), cmax(0.3, 0.9, 0.2);
    g.build(cmin, cmax, Vector3(-0.26, 0.04, -0.16), Vector3(0.26, 0.42, 0.16), 0.06);
    double dt = 1.2e-3;
    for (int s = 0; s < 400; s++) g.step(dt);
    double surf = 0; for (auto& p : g.pos) surf = std::max(surf, (double)p.y);
    FluidBall b; b.active = true; b.radius = 0.09; b.mass = 1.2;
    b.pos = Vector3(0, surf + 0.2, 0); b.vel = Vector3(0, -3.0, 0);
    g.balls.push_back(b);
    for (int s = 0; s < 700; s++) g.step(dt);
    return g.balls[0].pos.y;
}

int main() {
    // A) the generalised-Newtonian viscosity has the right shape
    Rheology ob = Rheology::oobleck(), wa = Rheology::water(), ke = Rheology::ketchup();
    CHECK(ob.viscosity(50) > ob.viscosity(1) * 5);        // shear-thickening: rises sharply
    CHECK_NEAR(wa.viscosity(50), wa.viscosity(1), 1e-6);  // Newtonian: flat
    CHECK(ke.viscosity(1) > ke.viscosity(50) * 5);        // shear-thinning: falls sharply
    CHECK(ke.viscosity(0.01) > 50);                       // yield stress: huge at low shear
    CHECK(wa.viscosity(0.0) >= wa.muMin);                 // regularised, finite at zero shear

    // B) a pool is stable and conservative: it settles, stays inside the box, no NaNs
    {
        SPHFluid g; Vector3 cmin(-0.3, 0, -0.2), cmax(0.3, 0.9, 0.2);
        g.build(cmin, cmax, Vector3(-0.26, 0.04, -0.16), Vector3(0.26, 0.5, 0.16), 0.06);
        int n0 = g.nParticles(); CHECK(n0 > 100); CHECK(g.rho0 > 0);
        for (int s = 0; s < 1200; s++) g.step(1.2e-3);
        CHECK(g.nParticles() == n0);                      // particle count conserved
        double vmax = 0; bool finite = true, inside = true;
        for (int i = 0; i < n0; i++) {
            double v = std::sqrt(g.vel[i].squareMagnitude());
            vmax = std::max(vmax, v);
            if (!std::isfinite(v) || !std::isfinite(g.rho[i])) finite = false;
            for (int a = 0; a < 3; a++) if (g.pos[i][a] < cmin[a] - 1e-6 || g.pos[i][a] > cmax[a] + 1e-6) inside = false;
        }
        CHECK(finite);
        CHECK(inside);
        CHECK(vmax < 1.0);                                // calmed down (was a dropped block)
    }

    // C) behavioural: a shear-thickening pool holds an impactor higher than water
    CHECK(dropRestY(false) > dropRestY(true));

    // D) runtime add/remove keeps arrays in sync; a cylinder excludes particles
    {
        SPHFluid g; Vector3 cmin(-0.3, 0, -0.3), cmax(0.3, 0.6, 0.3);
        g.build(cmin, cmax, Vector3(-0.26, 0.04, -0.26), Vector3(-0.14, 0.24, -0.14), 0.06);   // fill clear of centre
        int n0 = g.nParticles();
        g.addParticle(Vector3(0, 0.5, 0), Vector3());
        CHECK(g.nParticles() == n0 + 1);
        CHECK((int)g.rho.size() == n0 + 1 && (int)g.vel.size() == n0 + 1);
        g.removeParticle(0);
        CHECK(g.nParticles() == n0);
        // a particle driven into the cylinder footprint is pushed back onto its surface
        SPHCylinder c; c.cx = 0; c.cz = 0; c.radius = 0.12; c.y0 = 0; c.y1 = 1e9; g.cylinders.push_back(c);
        g.addParticle(Vector3(0.04, 0.2, 0.0), Vector3(-1.0, 0, 0));   // inside the footprint, moving inward
        g.step(1.2e-3);
        int last = g.nParticles() - 1;
        real dx = g.pos[last].x, dz = g.pos[last].z;
        CHECK(std::sqrt(dx * dx + dz * dz) >= c.radius - 1e-3);        // ejected to the wall, not left inside
    }

    return test::report("sph");
}
