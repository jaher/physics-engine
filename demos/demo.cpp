// A headless demo that exercises the engine across its layers and writes an SVG
// dashboard (four panels): a bouncing particle, projectiles with/without drag, a
// rigid box dropped onto the ground, and a rod pendulum. Run: ./demo out.svg
#include "phys/phys.h"
#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
using namespace phys;

struct Poly { std::vector<std::pair<double, double>> pts; std::string colour; std::string label; };

// Emit one panel: axes + coloured polylines, given data-space bounds.
static std::string panel(double px, double py, double pw, double ph, const std::string& title,
                         const std::string& xlab, const std::string& ylab,
                         double x0, double x1, double y0, double y1, const std::vector<Poly>& polys) {
    static int uid = 0; int id = ++uid;
    std::ostringstream s;
    auto X = [&](double x) { return px + (x - x0) / (x1 - x0) * pw; };
    auto Y = [&](double y) { return py + ph - (y - y0) / (y1 - y0) * ph; };
    s << "<g font-family='Georgia,serif'>";
    s << "<clipPath id='clip" << id << "'><rect x='" << px << "' y='" << py << "' width='" << pw << "' height='" << ph << "'/></clipPath>";
    s << "<rect x='" << px << "' y='" << py << "' width='" << pw << "' height='" << ph
      << "' fill='#fbf8f0' stroke='#3a342a' stroke-width='1.5'/>";
    // zero line if in range
    if (y0 < 0 && y1 > 0) s << "<line x1='" << px << "' y1='" << Y(0) << "' x2='" << (px + pw) << "' y2='" << Y(0)
                            << "' stroke='#c9bfa8' stroke-width='1'/>";
    s << "<g clip-path='url(#clip" << id << ")'>";
    for (const auto& p : polys) {
        s << "<polyline fill='none' stroke='" << p.colour << "' stroke-width='2' points='";
        for (auto& pt : p.pts) s << X(pt.first) << "," << Y(pt.second) << " ";
        s << "'/>";
    }
    s << "</g>";
    s << "<text x='" << (px + 8) << "' y='" << (py + 20) << "' font-size='16' fill='#2a2620'>" << title << "</text>";
    s << "<text x='" << (px + pw - 6) << "' y='" << (py + ph - 6) << "' font-size='11' fill='#6a6252' text-anchor='end'>" << xlab << "</text>";
    s << "<text x='" << (px + 6) << "' y='" << (py + ph - 8) << "' font-size='11' fill='#6a6252'>" << ylab << "</text>";
    double ly = py + 34;
    for (const auto& p : polys) {
        s << "<line x1='" << (px + 10) << "' y1='" << ly << "' x2='" << (px + 30) << "' y2='" << ly
          << "' stroke='" << p.colour << "' stroke-width='3'/>";
        s << "<text x='" << (px + 35) << "' y='" << (ly + 4) << "' font-size='12' fill='#2a2620'>" << p.label << "</text>";
        ly += 16;
    }
    s << "</g>";
    return s.str();
}

int main(int argc, char** argv) {
    const char* out = argc > 1 ? argv[1] : "physics_demo.svg";

    // 1) Bouncing particle (ground restitution 0.7).
    Poly ball; ball.colour = "#2f4b8c"; ball.label = "height (m)";
    { Particle p; p.setMass(1); p.setDamping(1); p.setAcceleration(0, -9.81, 0); p.setPosition(0, 5, 0);
      ParticleWorld w(1); w.getParticles().push_back(&p);
      GroundContacts g; g.restitution = 0.7; g.init(&w.getParticles()); w.getContactGenerators().push_back(&g);
      for (int i = 0; i < 1440; i++) { w.startFrame(); w.runPhysics(1.0 / 240); ball.pts.push_back({i / 240.0, p.getPosition().y}); } }

    // 2) Projectiles: ideal parabola vs quadratic drag.
    Poly ideal, dragged; ideal.colour = "#2f8c4b"; ideal.label = "no drag"; dragged.colour = "#b04024"; dragged.label = "with drag";
    for (int mode = 0; mode < 2; mode++) {
        Particle q; q.setMass(1); q.setDamping(1); q.setAcceleration(0, -9.81, 0); q.setPosition(0, 0, 0); q.setVelocity(22, 20, 0);
        ParticleDrag drag(0.15, 0.045);
        ParticleWorld w(1); w.getParticles().push_back(&q);
        if (mode == 1) w.getForceRegistry().add(&q, &drag);
        for (int i = 0; i < 2000; i++) { w.startFrame(); w.runPhysics(1.0 / 120);
            Vector3 pos = q.getPosition(); if (pos.y < 0) break;
            (mode == 0 ? ideal : dragged).pts.push_back({pos.x, pos.y}); }
    }

    // 3) Rigid box dropped onto the ground: height and kinetic energy vs time.
    Poly boxH, boxE; boxH.colour = "#2f4b8c"; boxH.label = "height (m)"; boxE.colour = "#b04024"; boxE.label = "kinetic energy (J)";
    { struct BoxPlaneGen : ContactGenerator {
          CollisionBox* box; CollisionPlane plane; real fr, re;
          unsigned addContact(Contact* c, unsigned limit) const override { box->calculateInternals();
              CollisionData d; d.contactArray = c; d.reset(limit); d.friction = fr; d.restitution = re;
              CollisionDetector::boxAndHalfSpace(*box, plane, &d); return d.contactCount; } };
      RigidBody b; b.setMass(1); Matrix3 I; I.setBlockInertiaTensor(Vector3(0.5, 0.5, 0.5), 1); b.setInertiaTensor(I);
      b.setPosition(0, 4, 0); b.setOrientation(Quaternion(0.96, 0.2, 0.15, 0)); b.setAcceleration(0, -9.81, 0);
      b.setDamping(1, 1); b.calculateDerivedData();
      CollisionBox box; box.body = &b; box.halfSize = Vector3(0.5, 0.5, 0.5);
      CollisionPlane ground; ground.direction = Vector3(0, 1, 0); ground.offset = 0;
      BoxPlaneGen gen; gen.box = &box; gen.plane = ground; gen.fr = 0.4; gen.re = 0.45;
      World w(20); w.getRigidBodies().push_back(&b); w.getContactGenerators().push_back(&gen);
      for (int i = 0; i < 1000; i++) { w.startFrame(); w.runPhysics(1.0 / 200);
          boxH.pts.push_back({i / 200.0, b.getPosition().y});
          boxE.pts.push_back({i / 200.0, b.getKineticEnergy()}); }
      // scale energy to share the height axis (impact spikes reach tens of joules)
      double maxE = 1e-9; for (auto& p : boxE.pts) maxE = std::max(maxE, p.second);
      for (auto& p : boxE.pts) p.second *= 4.0 / maxE;
      char lbl[64]; std::snprintf(lbl, sizeof(lbl), "kinetic energy (×%.2f J)", maxE / 4.0); boxE.label = lbl; }

    // 4) Rod pendulum: the bob's traced path (x,y).
    Poly path; path.colour = "#7a3fa0"; path.label = "bob path";
    { Particle p; p.setMass(1); p.setDamping(0.999); p.setAcceleration(0, -9.81, 0); p.setPosition(1.6, 0, 0);
      ParticleRodConstraint rod; rod.particle = &p; rod.anchor = Vector3(0, 0, 0); rod.length = 1.6;
      ParticleWorld w(1); w.getParticles().push_back(&p); w.getContactGenerators().push_back(&rod);
      for (int i = 0; i < 2400; i++) { w.startFrame(); w.runPhysics(1.0 / 240);
          path.pts.push_back({p.getPosition().x, p.getPosition().y}); } }

    std::ostringstream svg;
    svg << "<svg xmlns='http://www.w3.org/2000/svg' width='1000' height='740' viewBox='0 0 1000 740'>";
    svg << "<rect width='1000' height='740' fill='#efe7d3'/>";
    svg << "<text x='500' y='30' font-family='Georgia,serif' font-size='22' fill='#2a2620' text-anchor='middle'>"
        << "Millington physics engine — verification dashboard</text>";
    svg << panel(30, 50, 450, 300, "Bouncing particle (restitution 0.7)", "time (s)", "y", 0, 6, 0, 5.2, {ball});
    svg << panel(520, 50, 450, 300, "Projectile: ideal vs quadratic drag", "x (m)", "y", 0, 55, 0, 24, {ideal, dragged});
    svg << panel(30, 390, 450, 320, "Rigid box drop: height and energy", "time (s)", "", 0, 5, 0, 4.2, {boxH, boxE});
    svg << panel(520, 390, 450, 320, "Rod pendulum bob path", "x (m)", "y", -1.9, 1.9, -1.9, 0.6, {path});
    svg << "</svg>";

    FILE* f = fopen(out, "w");
    if (!f) { std::perror("fopen"); return 1; }
    std::string str = svg.str(); fwrite(str.data(), 1, str.size(), f); fclose(f);
    std::printf("wrote %s  (%zu bytes)\n", out, str.size());
    std::printf("bounce peaks decay; drag shortens range; box settles with energy→0; pendulum arcs.\n");
    return 0;
}
