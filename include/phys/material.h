// Surface materials and how two of them combine into the effective coefficients
// of a contact — the "ContactMaterial / combine-mode" model used by PhysX, Unity
// and Bullet. A Material carries the four Coulomb coefficients the resolver needs:
// sliding friction, restitution, rolling friction and spinning (torsional)
// friction. `combine(a, b, mode)` merges two surfaces into the per-contact values.
#pragma once
#include "precision.h"

namespace phys {

// How two per-body coefficients are merged into a single contact coefficient.
enum class CombineMode { Average, Min, Max, Multiply };

inline real combineValue(real a, real b, CombineMode mode) {
    switch (mode) {
        case CombineMode::Min:      return a < b ? a : b;
        case CombineMode::Max:      return a > b ? a : b;
        case CombineMode::Multiply: return a * b;
        case CombineMode::Average:  default: return (a + b) * (real)0.5;
    }
}

struct Material {
    real friction        = (real)0.5;   // tangential (sliding) Coulomb coefficient
    real restitution     = (real)0.0;   // bounciness in [0,1]
    real rollingFriction = (real)0.0;   // resists rolling (torque opposing roll)
    real spinFriction    = (real)0.0;   // resists spin about the contact normal

    Material() {}
    Material(real f, real r, real roll = 0, real spin = 0)
        : friction(f), restitution(r), rollingFriction(roll), spinFriction(spin) {}
};

// Effective per-contact coefficients produced by merging two surface materials.
// The same combine mode is applied to every coefficient.
inline Material combine(const Material& a, const Material& b, CombineMode mode = CombineMode::Average) {
    return Material(
        combineValue(a.friction,        b.friction,        mode),
        combineValue(a.restitution,     b.restitution,     mode),
        combineValue(a.rollingFriction, b.rollingFriction, mode),
        combineValue(a.spinFriction,    b.spinFriction,    mode));
}

} // namespace phys
