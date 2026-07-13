// Forward-mode automatic differentiation: a Dual number carries (value, derivative)
// and propagates the chain rule through every operator and elementary function, so
// a simulation *templated on its scalar type* differentiates itself EXACTLY — no
// finite differences, no tuning. This is the "differentiable physics" that makes
// Brax useful: run the sim in Dual, read the analytic gradient out of `.d`.
#pragma once
#include "precision.h"
#include <cmath>

namespace phys {

// A dual number x = v + d·ε with ε² = 0. `v` is the value f(x₀); `d` carries the
// derivative f'(x₀). A `real` promotes to a *constant* Dual (derivative 0), so
// scalars in a formula are handled automatically; seed the independent variable
// with `variable(x)` (derivative 1).
struct Dual {
    real v, d;
    Dual() : v(0), d(0) {}
    Dual(real value) : v(value), d(0) {}            // a constant: dc/dx = 0
    Dual(real value, real deriv) : v(value), d(deriv) {}
};

inline Dual variable(real x) { return Dual(x, 1); }  // the thing we differentiate w.r.t.

// ---- arithmetic: each rule is the derivative of the operation ------------------
inline Dual operator+(const Dual& a, const Dual& b) { return Dual(a.v + b.v, a.d + b.d); }
inline Dual operator-(const Dual& a, const Dual& b) { return Dual(a.v - b.v, a.d - b.d); }
inline Dual operator-(const Dual& a)                { return Dual(-a.v, -a.d); }
inline Dual operator*(const Dual& a, const Dual& b) { return Dual(a.v * b.v, a.d * b.v + a.v * b.d); }
inline Dual operator/(const Dual& a, const Dual& b) {           // quotient rule
    real inv = ((real)1) / b.v;
    return Dual(a.v * inv, (a.d * b.v - a.v * b.d) * inv * inv);
}
inline Dual& operator+=(Dual& a, const Dual& b) { a = a + b; return a; }
inline Dual& operator-=(Dual& a, const Dual& b) { a = a - b; return a; }
inline Dual& operator*=(Dual& a, const Dual& b) { a = a * b; return a; }
inline Dual& operator/=(Dual& a, const Dual& b) { a = a / b; return a; }

// comparisons act on the value so control flow tracks the primal simulation
inline bool operator< (const Dual& a, const Dual& b) { return a.v <  b.v; }
inline bool operator> (const Dual& a, const Dual& b) { return a.v >  b.v; }
inline bool operator<=(const Dual& a, const Dual& b) { return a.v <= b.v; }
inline bool operator>=(const Dual& a, const Dual& b) { return a.v >= b.v; }
inline bool operator==(const Dual& a, const Dual& b) { return a.v == b.v; }

// ---- elementary functions: value, then value·(local derivative) ----------------
// The `using std::…` lets a templated simulation call sin/cos/… unqualified and
// get std::sin for `real` and phys::sin for `Dual` from one overload set.
using std::sin; using std::cos; using std::exp; using std::sqrt; using std::pow; using std::fabs;
inline Dual sin (const Dual& a) { return Dual(std::sin(a.v),  std::cos(a.v) * a.d); }
inline Dual cos (const Dual& a) { return Dual(std::cos(a.v), -std::sin(a.v) * a.d); }
inline Dual exp (const Dual& a) { real e = std::exp(a.v);  return Dual(e, e * a.d); }
inline Dual sqrt(const Dual& a) { real s = std::sqrt(a.v); return Dual(s, a.d / (2 * s)); }
inline Dual pow (const Dual& a, real p) { return Dual(std::pow(a.v, p), p * std::pow(a.v, p - 1) * a.d); }
inline Dual abs (const Dual& a) { return Dual(std::fabs(a.v), (a.v < 0 ? -a.d : a.d)); }
inline Dual fabs(const Dual& a) { return abs(a); }

// ================ tiny simulations, templated on the scalar type ================
// Instantiate with `real` for the ordinary result, or with `Dual` (seeding one
// input via `variable`) to read the analytic gradient straight out of `.d`.

// Ideal projectile launched at `angle` with speed v0 under gravity g → horizontal
// range. Differentiating in Dual gives the exact d(range)/d(angle).
template<class T>
inline T projectileRange(T angle, real v0, real g) {
    T vx = v0 * cos(angle);              // horizontal speed
    T vy = v0 * sin(angle);              // launch vertical speed
    T timeOfFlight = (2 * vy) * (1 / g); // up-and-back under gravity
    return vx * timeOfFlight;            // range = vx · t
}

// A damped spring released from x0, advanced `steps` of semi-implicit Euler.
// Returns the final displacement as a function of stiffness `k`; templating on the
// scalar type differentiates the *discrete* simulation w.r.t. stiffness — the
// gradient flows through every one of the time steps.
template<class T>
inline T springSettle(T k, real m, real c, real x0, int steps, real dt) {
    T x = T(x0), v = T((real)0);
    for (int s = 0; s < steps; s++) {
        T a = (-(k * x) - c * v) * (1 / m);   // F = -k·x - c·v ,  a = F/m
        v = v + a * dt;
        x = x + v * dt;
    }
    return x;
}

} // namespace phys
