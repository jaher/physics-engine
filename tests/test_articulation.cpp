// Featherstone reduced-coordinate articulated-body dynamics (ABA). We check it
// against closed-form results: the compound-pendulum period, the released-link
// angular acceleration m g l / I_pivot, energy conservation of a free 2-link
// chain, and a gravity-cancelling holding torque.
#include "phys/articulation.h"
#include "check.h"
#include <cmath>
using namespace phys;

static const real G = 9.81;

// A thin rod of length L, mass m, lying along +x with its inertia about the COM
// (I_zz = m L²/12). `comAlong` places the COM offset from the joint anchor.
static Matrix3 rodInertia(real m, real L) {
    real It = m * L * L / 12;
    Matrix3 I; I.setInertiaTensorCoeffs(1e-4 * m, It, It);      // tiny Ixx avoids a zero row
    return I;
}

int main() {
    // ---------------------------------------------------------------------- (a)
    // Single revolute link as a compound pendulum hanging straight down: measure
    // the small-oscillation period and compare to T = 2π√(I_pivot/(m g l)).
    {
        real m = 1.0, L = 1.0, l = 0.5;                        // COM at rod centre
        real Ipivot = m * L * L / 12 + m * l * l;              // = m L²/3
        real Tref = 2 * real_pi * std::sqrt(Ipivot / (m * G * l));

        Articulation art;
        art.setGravity(Vector3(0, -G, 0));
        int link = art.addLink(-1, JointType::REVOLUTE, Vector3(0, 0, 1),
                               m, rodInertia(m, L), Vector3(0, -l, 0), Vector3(0, 0, 0));
        art.setJointState(link, 0.05, 0);                      // small release angle

        // Integrate and time a full period from two successive upward zero-crossings.
        real dt = 2e-4, prev = art.q(link), t = 0;
        real firstCross = -1, secondCross = -1; real prevT = 0;
        for (int s = 0; s < 200000 && secondCross < 0; s++) {
            art.step(dt); t += dt;
            real cur = art.q(link);
            if (prev < 0 && cur >= 0) {                        // rising zero-crossing
                real tc = t - dt * cur / (cur - prev);         // linear interp
                if (firstCross < 0) firstCross = tc; else secondCross = tc;
            }
            prev = cur; prevT = t;
        }
        (void)prevT;
        CHECK(secondCross > 0);
        real Tmeas = secondCross - firstCross;
        CHECK_NEAR(Tmeas, Tref, 0.03 * Tref);                  // within 3%
    }

    // ---------------------------------------------------------------------- (b)
    // Torque-free horizontal link released from rest: initial angular accel must
    // be −m g l / I_pivot (COM at +x falls → negative rotation about +z).
    {
        real m = 2.0, L = 1.4, l = 0.7;
        real Ipivot = m * L * L / 12 + m * l * l;
        Articulation art;
        art.setGravity(Vector3(0, -G, 0));
        int link = art.addLink(-1, JointType::REVOLUTE, Vector3(0, 0, 1),
                               m, rodInertia(m, L), Vector3(l, 0, 0), Vector3(0, 0, 0));
        art.step(1e-6);                                        // one solve from rest
        CHECK_NEAR(art.qddot(link), -m * G * l / Ipivot, 1e-6);
    }

    // ---------------------------------------------------------------------- (c)
    // Free (zero-torque) 2-link chain run ~2 s conserves total energy to a few %.
    {
        real m = 1.0, L = 1.0, l = 0.5;
        Articulation art;
        art.setGravity(Vector3(0, -G, 0));
        int a = art.addLink(-1, JointType::REVOLUTE, Vector3(0, 0, 1),
                            m, rodInertia(m, L), Vector3(l, 0, 0), Vector3(0, 0, 0));
        int b = art.addLink(a, JointType::REVOLUTE, Vector3(0, 0, 1),
                            m, rodInertia(m, L), Vector3(l, 0, 0), Vector3(L, 0, 0));
        art.setJointState(a, 0.8, 0);                          // start from a bent, raised pose
        art.setJointState(b, -0.6, 0);
        (void)b;

        real E0 = art.totalEnergy();
        real dt = 1e-4;
        for (int s = 0; s < 20000; s++) art.step(dt);          // 2 s
        real E1 = art.totalEnergy();
        CHECK_NEAR(E1, E0, 0.03 * std::fabs(E0));              // < 3% drift
    }

    // ---------------------------------------------------------------------- (d)
    // A holding torque τ = m g l exactly cancels gravity on a horizontal link,
    // which should therefore stay put.
    {
        real m = 1.5, L = 1.0, l = 0.5;
        Articulation art;
        art.setGravity(Vector3(0, -G, 0));
        int link = art.addLink(-1, JointType::REVOLUTE, Vector3(0, 0, 1),
                               m, rodInertia(m, L), Vector3(l, 0, 0), Vector3(0, 0, 0));
        art.setJointTorque(link, m * G * l);                   // balances weight at q=0
        real dt = 5e-4;
        for (int s = 0; s < 1000; s++) art.step(dt);           // 0.5 s
        CHECK_NEAR(art.q(link), 0.0, 1e-3);
        CHECK_NEAR(art.qd(link), 0.0, 1e-3);
    }

    return test::report("articulation");
}
