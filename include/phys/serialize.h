// Snapshot & restore of a set of rigid bodies to/from a flat byte buffer — the
// deterministic-rollback / networking primitive engines expose (PhysX serialization,
// MuJoCo mjData save/load, Bullet world serializer). Each body contributes its
// position, orientation quaternion, and linear + angular velocity. The doubles are
// copied verbatim, so position and both velocities round-trip *bitwise*; the
// orientation is written back through the body's normalising setter (it stays a
// unit quaternion, reproduced to full double precision up to renormalisation).
#pragma once
#include "body.h"
#include <vector>
#include <cstring>
#include <cstdint>

namespace phys {

struct Snapshot {
    // 13 doubles per body: pos(3) + quat(4) + linVel(3) + angVel(3)
    static const int kDoublesPerBody = 13;
    static const uint32_t kMagic = 0x50485331;   // "PHS1"

    // Serialise the bodies' dynamic state into a byte buffer.
    static std::vector<unsigned char> snapshot(const std::vector<RigidBody*>& bodies) {
        uint32_t count = (uint32_t)bodies.size();
        std::vector<unsigned char> buf(8 + (size_t)count * kDoublesPerBody * sizeof(double));
        std::memcpy(buf.data() + 0, &kMagic, 4);
        std::memcpy(buf.data() + 4, &count, 4);
        double* d = reinterpret_cast<double*>(buf.data() + 8);
        for (uint32_t b = 0; b < count; b++) {
            const RigidBody* body = bodies[b];
            Vector3 p = body->getPosition();
            Quaternion q = body->getOrientation();
            Vector3 v = body->getVelocity();
            Vector3 w = body->getRotation();
            *d++ = p.x; *d++ = p.y; *d++ = p.z;
            *d++ = q.r; *d++ = q.i; *d++ = q.j; *d++ = q.k;
            *d++ = v.x; *d++ = v.y; *d++ = v.z;
            *d++ = w.x; *d++ = w.y; *d++ = w.z;
        }
        return buf;
    }

    // Restore state previously produced by snapshot() onto the *same-sized* body
    // list (same order). Returns false on a malformed buffer or a size mismatch.
    static bool restore(const std::vector<unsigned char>& buf, const std::vector<RigidBody*>& bodies) {
        if (buf.size() < 8) return false;
        uint32_t magic = 0, count = 0;
        std::memcpy(&magic, buf.data() + 0, 4);
        std::memcpy(&count, buf.data() + 4, 4);
        if (magic != kMagic) return false;
        if (count != bodies.size()) return false;
        if (buf.size() != 8 + (size_t)count * kDoublesPerBody * sizeof(double)) return false;
        const double* d = reinterpret_cast<const double*>(buf.data() + 8);
        for (uint32_t b = 0; b < count; b++) {
            RigidBody* body = bodies[b];
            Vector3 p(d[0], d[1], d[2]);
            Quaternion q(d[3], d[4], d[5], d[6]);
            Vector3 v(d[7], d[8], d[9]);
            Vector3 w(d[10], d[11], d[12]);
            d += kDoublesPerBody;
            body->setPosition(p);
            body->setOrientation(q);       // normalises (keeps it a unit quaternion)
            body->setVelocity(v);
            body->setRotation(w);
            body->calculateDerivedData();  // rebuild transform + world inverse inertia
        }
        return true;
    }

    // Convenience: number of bytes a snapshot of `n` bodies occupies.
    static size_t byteSize(size_t n) { return 8 + n * kDoublesPerBody * sizeof(double); }
};

} // namespace phys
