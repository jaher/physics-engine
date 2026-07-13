// A small deterministic random source for demos (Millington appendix) — vectors,
// scalars and quaternions, seeded reproducibly.
#pragma once
#include "core.h"
#include <cstdint>

namespace phys {

class Random {
    uint32_t s;
public:
    Random(uint32_t seed = 0x1234abcd) { this->seed(seed); }
    void seed(uint32_t s_) { s = s_ ? s_ : 1; }
    uint32_t randomBits() { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return s; }   // xorshift32
    real randomReal() { return (randomBits() & 0xffffff) / (real)0x1000000; }        // [0,1)
    real randomReal(real min, real max) { return min + randomReal() * (max - min); }
    real randomReal(real scale) { return randomReal() * scale; }
    unsigned randomInt(unsigned max) { return randomBits() % max; }
    Vector3 randomVector(real scale) {
        return Vector3(randomReal(-scale, scale), randomReal(-scale, scale), randomReal(-scale, scale));
    }
    Vector3 randomVector(const Vector3& min, const Vector3& max) {
        return Vector3(randomReal(min.x, max.x), randomReal(min.y, max.y), randomReal(min.z, max.z));
    }
    Quaternion randomQuaternion() {
        Quaternion q(randomReal(), randomReal(), randomReal(), randomReal());
        q.normalise(); return q;
    }
};

} // namespace phys
