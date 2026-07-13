// Coupled spring-mass lattice → normal modes and standing-wave harmonics.
// A line of point masses (phys::Particle) joined by Hooke springs
// (phys::ParticleSpring, one per direction so Newton's third law holds) with the
// two end nodes held fixed. Excited in a single normal mode the chain rings as a
// standing wave; its mode frequencies follow the analytic monatomic-lattice
// dispersion relation
//     ω_n = 2·sqrt(k/m)·sin( nπ / (2(P+1)) ),   n = 1 … P,
// with interior mode shape  a_j = sin( nπ j /(P+1) ),  j = 1 … P.
// Longitudinal motion is exactly linear (the rest length cancels); transverse
// motion needs the springs pre-tensioned (restLength < spacing) and is linear to
// first order, giving the familiar harmonic series of a plucked string.
#pragma once
#include "particle.h"
#include "pfgen.h"
#include <vector>

namespace phys {

class SpringChain {
public:
    std::vector<Particle> node;          // node[0] and node[P+1] are the fixed walls
    std::vector<ParticleSpring> springs; // two per link (one force per endpoint)
    ParticleForceRegistry registry;

    Vector3 origin, axis = Vector3(1, 0, 0);
    real spacing = 0, restLength = 0, k = 0, m = 1;
    int P = 0;                           // number of movable interior masses

    SpringChain() = default;
    SpringChain(const SpringChain&) = delete;            // registry holds interior pointers
    SpringChain& operator=(const SpringChain&) = delete;

    // Build P interior masses between fixed ends 0 and P+1, evenly spaced along
    // `axis`. restLength < spacing pre-tensions the springs (needed for crisp
    // transverse harmonics); restLength == spacing leaves them slack (fine for
    // longitudinal modes, which don't depend on tension).
    void build(int P_, real spacing_, real k_, real m_, const Vector3& origin_,
               const Vector3& axis_ = Vector3(1, 0, 0), real restLength_ = -1, real damping = 1.0) {
        P = P_; spacing = spacing_; k = k_; m = m_; origin = origin_; axis = axis_.unit();
        restLength = restLength_ < 0 ? spacing_ : restLength_;
        node.assign(P + 2, Particle());
        for (int i = 0; i <= P + 1; i++) {
            node[i].setPosition(origin + axis * (spacing * i));
            node[i].setVelocity(Vector3());
            node[i].setAcceleration(Vector3());        // no gravity: pure spring dynamics
            node[i].setDamping(damping);
            if (i == 0 || i == P + 1) node[i].setInverseMass(0);   // immovable walls
            else node[i].setMass(m);
        }
        springs.clear(); springs.reserve(2 * (P + 1));             // reserve → stable addresses
        for (int i = 0; i <= P; i++) {
            springs.emplace_back(&node[i + 1], k, restLength);     // force on node[i]
            springs.emplace_back(&node[i], k, restLength);         // force on node[i+1]
        }
        registry.clear();
        for (int i = 0; i <= P; i++) {
            registry.add(&node[i], &springs[2 * i]);
            registry.add(&node[i + 1], &springs[2 * i + 1]);
        }
    }

    // One symplectic (semi-implicit) Euler step of the whole chain: accumulate the
    // Hooke forces (ParticleSpring) at the current positions, advance each velocity,
    // then advance each position with the *updated* velocity. This ordering is
    // symplectic, so a ringing chain keeps its energy bounded instead of pumping it
    // up the way plain explicit Euler would.
    void step(real dt) {
        for (auto& p : node) p.clearAccumulator();
        registry.updateForces(dt);
        for (auto& p : node) {
            if (!p.hasFiniteMass()) continue;
            Vector3 acc = p.getAcceleration() + p.getForceAccum() * p.getInverseMass();
            Vector3 v = p.getVelocity() + acc * dt;
            v *= real_pow(p.getDamping(), dt);
            p.setVelocity(v);
            p.setPosition(p.getPosition() + v * dt);
        }
    }

    // --- analytic normal-mode helpers ---
    Vector3 equilibrium(int j) const { return origin + axis * (spacing * j); }
    real modeShape(int n, int j) const { return real_sin(n * real_pi * j / (P + 1)); }
    // effective coupling: k for longitudinal, tension/spacing for transverse
    real modeOmega(int n, bool transverse = false) const {
        real keff = transverse ? k * (spacing - restLength) / spacing : k;
        return 2 * real_sqrt(keff / m) * real_sin(n * real_pi / (2 * (P + 1)));
    }
    real modePeriod(int n, bool transverse = false) const { return 2 * real_pi / modeOmega(n, transverse); }

    // Set the chain to pure mode n (zero velocity): displacement along `dir`,
    // amplitude `amp`. dir == axis gives a longitudinal mode; a perpendicular dir
    // gives a transverse standing wave.
    void setMode(int n, real amp, const Vector3& dir) {
        Vector3 d = dir.unit();
        for (int j = 1; j <= P; j++) {
            node[j].setPosition(equilibrium(j) + d * (amp * modeShape(n, j)));
            node[j].setVelocity(Vector3());
        }
    }

    // Displacement of interior node j from its equilibrium.
    Vector3 displacement(int j) const { return node[j].getPosition() - equilibrium(j); }

    // Total mechanical energy (kinetic + spring potential) — conserved when damping = 1.
    real energy() const {
        real E = 0;
        for (const auto& p : node) E += p.getKineticEnergy();
        for (int i = 0; i <= P; i++) {
            real x = (node[i + 1].getPosition() - node[i].getPosition()).magnitude() - restLength;
            E += (real)0.5 * k * x * x;
        }
        return E;
    }
};

} // namespace phys
