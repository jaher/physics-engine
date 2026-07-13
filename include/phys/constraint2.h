// Soft constraints and actuation the core joint zoo lacks: a CFM/ERP soft
// distance link (springy), a breakable joint, a conveyor-belt surface, wind and
// radial force fields, rigid-body hydrodynamics (buoyancy + added mass + drag +
// lift), and an unconditionally-stable backward-Euler mass-spring integrator.
// Soft constraints follow the ERP/CFM (Baumgarte + constraint-force-mixing)
// formulation used by ODE/Bullet; the implicit spring is standard backward Euler.
#pragma once
#include "body.h"
#include <vector>
#include <functional>

namespace phys {

// --- shared soft-constraint helpers -------------------------------------------
// World-space velocity of a body point (linear + w×r).
inline Vector3 c2_pointVel(const RigidBody* b, const Vector3& worldPt) {
    if (!b) return Vector3();
    return b->getVelocity() + b->getRotation() % (worldPt - b->getPosition());
}
// Apply a linear impulse at a world point (updates linear + angular velocity).
inline void c2_applyImpulse(RigidBody* b, const Vector3& worldPt, const Vector3& impulse) {
    if (!b || b->getInverseMass() == 0) return;
    b->addVelocity(impulse * b->getInverseMass());
    Matrix3 iit; b->getInverseInertiaTensorWorld(&iit);
    b->addRotation(iit.transform((worldPt - b->getPosition()) % impulse));
    b->setAwake(true);
}
// Inverse effective mass along a unit direction n at the two anchors (linear +
// angular lever arms) — the K in λ = -(Jv+bias)/(K+CFM).
inline real c2_effInvMass(const RigidBody* b0, const Vector3& a0,
                          const RigidBody* b1, const Vector3& a1, const Vector3& n) {
    const RigidBody* bs[2] = {b0, b1}; Vector3 as[2] = {a0, a1};
    real k = 0;
    for (int i = 0; i < 2; i++) {
        const RigidBody* b = bs[i]; if (!b || b->getInverseMass() == 0) continue;
        Vector3 r = as[i] - b->getPosition();
        Matrix3 iit; b->getInverseInertiaTensorWorld(&iit);
        k += b->getInverseMass() + ((iit.transform(r % n)) % r) * n;   // ((I⁻¹(r×n))×r)·n
    }
    return k;
}

// --- 1. Soft distance constraint (CFM + ERP) ----------------------------------
// Holds |anchorA − anchorB| = restLength as a soft impulse each step. `softness`
// is the constraint-force-mixing term: 0 → a rigid rod, larger → a springier
// link that carries less force for the same error. `erp` bleeds position error
// out over ~1/erp steps. body[1]==nullptr anchors to a fixed world point.
class SoftDistanceConstraint {
public:
    RigidBody* body[2] = {nullptr, nullptr};
    Vector3 local[2];                       // anchor in each body's space (world if body null)
    real restLength = 1;
    real softness = 0;                      // CFM: 0 rigid, >0 springy (smaller restoring force)
    real erp = (real)0.2;                   // Baumgarte error-reduction fraction per step
    real lastImpulse = 0;                   // signed normal impulse from the last solve
    real lastForce = 0;                     // |impulse|/dt — the force the link had to carry

    void set(RigidBody* a, const Vector3& la, RigidBody* b, const Vector3& lb, real len) {
        body[0] = a; body[1] = b; local[0] = la; local[1] = lb; restLength = len;
    }
    Vector3 worldAnchor(int i) const {
        return body[i] ? body[i]->getPointInWorldSpace(local[i]) : local[i];
    }
    void solve(real dt) {
        Vector3 wa = worldAnchor(0), wb = worldAnchor(1);
        Vector3 d = wb - wa; real len = d.magnitude();
        if (len < (real)1e-9) { lastImpulse = lastForce = 0; return; }
        Vector3 n = d * (((real)1) / len);
        real C = len - restLength;                                     // position error
        real jv = n * (c2_pointVel(body[1], wb) - c2_pointVel(body[0], wa));
        real denom = c2_effInvMass(body[0], wa, body[1], wb, n) + softness;
        if (denom <= 0) { lastImpulse = lastForce = 0; return; }
        real bias = (erp / dt) * C;
        real lambda = -(jv + bias) / denom;
        lastImpulse = lambda; lastForce = real_abs(lambda) / dt;
        c2_applyImpulse(body[1], wb, n * lambda);
        c2_applyImpulse(body[0], wa, n * (-lambda));
    }
};

// --- 2. Breakable joint --------------------------------------------------------
// Wraps a soft constraint; the first time the force it must carry exceeds
// breakForce it detaches (stops solving) and fires onBreak exactly once.
class BreakableJoint {
public:
    SoftDistanceConstraint* constraint = nullptr;
    real breakForce = REAL_MAX;
    bool broken = false;
    std::function<void()> onBreak;          // fired once, when load first crosses breakForce

    void solve(real dt) {
        if (broken || !constraint) return;
        constraint->solve(dt);
        if (constraint->lastForce > breakForce) {
            broken = true;
            if (onBreak) onBreak();
        }
    }
    real load() const { return constraint ? constraint->lastForce : 0; }
};

// --- 3. Conveyor surface -------------------------------------------------------
// A belt: bodies resting on the plane (normal·com − offset in [0,maxHeight]) are
// dragged toward beltVelocity by a Coulomb-limited tangential impulse, so a box
// at rest spins up to belt speed and then rides along with it.
class ConveyorSurface {
public:
    Vector3 normal = Vector3(0, 1, 0);      // belt plane normal (world)
    real offset = 0;                        // plane: normal·x = offset
    Vector3 beltVelocity;                   // tangential surface velocity (world)
    real maxHeight = 1;                      // com within this band above the plane counts as resting
    real mu = 1;                            // friction coefficient limiting the drive
    Vector3 g = Vector3(0, (real)-9.81, 0); // gravity, for the friction normal-load estimate

    bool apply(RigidBody* body, real dt) {
        if (!body || !body->hasFiniteMass()) return false;
        real h = normal * body->getPosition() - offset;
        if (h < 0 || h > maxHeight) return false;                      // not on the belt
        Vector3 v = body->getVelocity();
        Vector3 vt = v - normal * (v * normal);                        // body tangential velocity
        Vector3 bt = beltVelocity - normal * (beltVelocity * normal);  // belt tangential velocity
        Vector3 slip = bt - vt; real slipMag = slip.magnitude();
        if (slipMag < (real)1e-9) return true;
        real maxDv = mu * real_abs(g * normal) * dt;                   // Coulomb cap on Δv this step
        real dv = slipMag < maxDv ? slipMag : maxDv;
        body->addVelocity(slip * (dv / slipMag));
        body->setAwake(true);
        return true;
    }
};

// --- 4a. Wind field ------------------------------------------------------------
// Uniform velocity field plus a sinusoidal gust; produces a linear+quadratic
// drag force pulling a body toward the local wind velocity (i.e. downwind).
class WindField {
public:
    Vector3 wind;                           // uniform base velocity
    Vector3 gustAxis = Vector3(1, 0, 0);
    real gustAmp = 0, gustFreq = 1;         // additive gust: gustAxis·gustAmp·sin(gustFreq·t)
    real linearDrag = (real)0.5;            // ∝ (wind − vel)
    real quadraticDrag = (real)0.1;         // ∝ |wind − vel|·(wind − vel)
    real time = 0;

    Vector3 velocityAt(real t) const { return wind + gustAxis * (gustAmp * real_sin(gustFreq * t)); }
    Vector3 dragForce(const Vector3& bodyVel, real t) const {
        Vector3 rel = velocityAt(t) - bodyVel;
        return rel * linearDrag + rel * (quadraticDrag * rel.magnitude());
    }
    void apply(RigidBody* body, real dt) {
        if (!body || !body->hasFiniteMass()) return;
        time += dt;
        body->addForce(dragForce(body->getVelocity(), time));
    }
};

// --- 4b. Radial field ----------------------------------------------------------
// Central attractor (strength>0) or repeller (strength<0) with either 1/r²
// (gravity-like) or linear (spring-like) falloff, clamped near the centre.
class RadialField {
public:
    Vector3 centre;
    real strength = 1;                      // >0 attracts toward centre, <0 repels
    bool inverseSquare = true;              // true: 1/r²   false: linear (∝ r)
    real minDist = (real)0.1;               // clamp radius to avoid the singularity
    real maxDist = REAL_MAX;                // beyond this the field is inert

    Vector3 forceAt(const Vector3& pos) const {
        Vector3 d = centre - pos; real r = d.magnitude();
        if (r < (real)1e-9 || r > maxDist) return Vector3();
        Vector3 dir = d * (((real)1) / r);                             // unit toward the centre
        real rc = r < minDist ? minDist : r;
        real mag = inverseSquare ? strength / (rc * rc) : strength * rc;
        return dir * mag;
    }
    void apply(RigidBody* body, real dt) {
        (void)dt;                               // stateless field; dt kept for a uniform apply() signature
        if (!body || !body->hasFiniteMass()) return;
        body->addForce(forceAt(body->getPosition()));
    }
};

// --- 5. Rigid-body hydrodynamics ----------------------------------------------
// Buoyancy + added mass + quadratic drag + a simple velocity-perpendicular lift
// for a body submerged in a fluid. A body denser than the fluid sinks to a
// bounded terminal velocity (weight − buoyancy balanced by quadratic drag).
class RigidHydrodynamics {
public:
    real fluidDensity = 1000;               // kg/m³
    real waterHeight = 0;                   // world y of the free surface
    real volume = (real)0.001;              // displaced volume when fully submerged
    real bodyHeight = (real)0.1;            // vertical extent, for partial submersion
    real dragCoeff = 1;                     // lumped quadratic-drag coefficient (½·Cd·A)
    real addedMass = (real)0.5;             // added-mass coefficient (× ρ·V_sub)
    real liftCoeff = 0;                     // simple lift coefficient (⟂ to velocity)
    Vector3 g = Vector3(0, (real)-9.81, 0);
    Vector3 prevVel; bool havePrev = false; // for the added-mass acceleration estimate

    real submergedFraction(real comY) const {
        real top = comY + bodyHeight * (real)0.5, bot = comY - bodyHeight * (real)0.5;
        if (bot >= waterHeight) return 0;                              // fully above water
        if (top <= waterHeight) return 1;                             // fully submerged
        return (waterHeight - bot) / bodyHeight;                       // partially
    }
    // Buoyancy + quadratic drag + lift (no added mass — that needs history/dt).
    Vector3 staticForce(const RigidBody* body) const {
        real frac = submergedFraction(body->getPosition().y);
        if (frac <= 0) return Vector3();
        real gmag = g.magnitude();
        Vector3 up = gmag > 0 ? g * ((real)-1 / gmag) : Vector3(0, 1, 0);
        real subVol = volume * frac;
        Vector3 f = up * (fluidDensity * subVol * gmag);              // buoyancy
        Vector3 v = body->getVelocity(); real s = v.magnitude();
        if (s > (real)1e-9) {
            Vector3 vh = v * (((real)1) / s);
            f += v * (-(real)0.5 * fluidDensity * dragCoeff * s * frac);   // quadratic drag
            if (liftCoeff != 0) {                                          // lift ⟂ v, in the v–up plane
                Vector3 perp = up - vh * (up * vh); real pl = perp.magnitude();
                if (pl > (real)1e-9)
                    f += perp * ((((real)1) / pl) * (real)0.5 * fluidDensity * liftCoeff * s * s * frac);
            }
        }
        return f;
    }
    void apply(RigidBody* body, real dt) {
        if (!body || !body->hasFiniteMass()) return;
        real frac = submergedFraction(body->getPosition().y);
        if (frac <= 0) { havePrev = false; return; }
        Vector3 f = staticForce(body);
        Vector3 v = body->getVelocity();
        if (havePrev && dt > 0) {                                     // added mass opposes acceleration
            Vector3 accel = (v - prevVel) * (((real)1) / dt);
            f += accel * (-(addedMass * fluidDensity * volume * frac));
        }
        prevVel = v; havePrev = true;
        body->addForce(f);
    }
};

// --- 6. Implicit (backward-Euler) mass-spring ---------------------------------
// Stiff 1-DOF mass-spring-damper about `anchor`:  m x'' = -k(x-anchor) - c x'.
// Solved implicitly for the new velocity, which is unconditionally stable — the
// energy stays bounded (in fact decays) for any dt, even huge k·dt²:
//   v⁺ = (m v - dt·k·(x-anchor)) / (m + dt·c + dt²·k),   x⁺ = x + dt·v⁺.
class ImplicitSpring1D {
public:
    real mass = 1, k = 1, c = 0;
    real x = 0, v = 0, anchor = 0;
    void step(real dt) {
        real xr = x - anchor;
        real denom = mass + dt * c + dt * dt * k;                     // always > 0
        v = (mass * v - dt * k * xr) / denom;
        x += dt * v;
    }
    real energy() const { real xr = x - anchor; return (real)0.5 * mass * v * v + (real)0.5 * k * xr * xr; }
};

// Explicit (forward) Euler on the same system — the unstable reference: for an
// undamped oscillator its energy grows every step, so a large k·dt² blows up.
class ExplicitSpring1D {
public:
    real mass = 1, k = 1, c = 0;
    real x = 0, v = 0, anchor = 0;
    void step(real dt) {
        real a = (-k * (x - anchor) - c * v) / mass;
        x += dt * v;                                                  // uses the old velocity
        v += dt * a;
    }
    real energy() const { real xr = x - anchor; return (real)0.5 * mass * v * v + (real)0.5 * k * xr * xr; }
};

// Backward Euler for a small chain of equal masses on a line with fixed ends,
// stiffness k between neighbours (and to each wall) and damping c. The implicit
// velocity update is a symmetric tridiagonal system solved with the Thomas
// algorithm each step, so the whole chain stays stable at large stiff dt.
class ImplicitSpringChain {
public:
    real mass = 1, k = 1, c = 0;
    std::vector<real> u, v;                 // offset from equilibrium, velocity per node
    void resize(int n) { u.assign(n, 0); v.assign(n, 0); }
    real energy() const {
        int n = (int)u.size(); real e = 0;
        for (int i = 0; i < n; i++) e += (real)0.5 * mass * v[i] * v[i];
        for (int i = 0; i <= n; i++) {      // n+1 springs (both ends anchored to walls at 0)
            real a = i > 0 ? u[i - 1] : 0, b = i < n ? u[i] : 0;
            e += (real)0.5 * k * (b - a) * (b - a);
        }
        return e;
    }
    void step(real dt) {
        int n = (int)u.size(); if (n == 0) return;
        real diag = mass + dt * c + 2 * dt * dt * k;                  // A_ii
        real off = -dt * dt * k;                                      // A_i,i±1
        std::vector<real> rhs(n), cp(n), dp(n);
        for (int i = 0; i < n; i++) {
            real ul = i > 0 ? u[i - 1] : 0, ur = i + 1 < n ? u[i + 1] : 0;
            rhs[i] = mass * v[i] - dt * k * (2 * u[i] - ul - ur);
        }
        cp[0] = off / diag; dp[0] = rhs[0] / diag;                    // Thomas forward sweep
        for (int i = 1; i < n; i++) {
            real m = diag - off * cp[i - 1];
            cp[i] = off / m;
            dp[i] = (rhs[i] - off * dp[i - 1]) / m;
        }
        v[n - 1] = dp[n - 1];                                         // back substitution
        for (int i = n - 2; i >= 0; i--) v[i] = dp[i] - cp[i] * v[i + 1];
        for (int i = 0; i < n; i++) u[i] += dt * v[i];
    }
};

} // namespace phys
