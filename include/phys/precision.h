// Precision configuration for the engine. Following Millington's "Game Physics
// Engine Development", all real arithmetic goes through `real` and the real_*
// wrappers so single/double precision is a one-line switch.
#pragma once
#include <cmath>
#include <cfloat>

namespace phys {

typedef double real;
#define REAL_MAX      DBL_MAX
#define REAL_EPSILON  DBL_EPSILON

const real real_pi = 3.14159265358979323846;

inline real real_sqrt(real x) { return std::sqrt(x); }
inline real real_abs (real x) { return std::fabs(x); }
inline real real_sin (real x) { return std::sin(x); }
inline real real_cos (real x) { return std::cos(x); }
inline real real_exp (real x) { return std::exp(x); }
inline real real_pow (real a, real b) { return std::pow(a, b); }
inline real real_fmod(real a, real b) { return std::fmod(a, b); }

} // namespace phys
