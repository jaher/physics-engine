// Heightfield terrain (ODE/PhysX/Bullet/Gazebo/Unity-style): a grid of height
// samples with bilinear interpolation, generating contacts for spheres, capsules
// and boxes, plus a terrain raycast.
#pragma once
#include "shapes.h"
#include <vector>

namespace phys {

class HeightField {
public:
    int nx = 0, nz = 0; real spacing = 1;
    Vector3 origin;                                        // corner at (0,0)
    std::vector<real> h;                                   // nx*nz samples
    void init(int nx_, int nz_, real spacing_, const Vector3& origin_) {
        nx = nx_; nz = nz_; spacing = spacing_; origin = origin_; h.assign(nx * nz, 0);
    }
    real& at(int x, int z) { return h[z * nx + x]; }
    // Press the surface down around world (wx,wz): an elliptical depression of
    // `depth` with radii (rx,rz), smooth falloff. Never deepens an existing
    // deeper print (max of depressions). This is how deformable snow works.
    void deform(real wx, real wz, real rx, real rz, real depth) {
        int x0 = (int)std::floor((wx - rx - origin.x) / spacing), x1 = (int)std::ceil((wx + rx - origin.x) / spacing);
        int z0 = (int)std::floor((wz - rz - origin.z) / spacing), z1 = (int)std::ceil((wz + rz - origin.z) / spacing);
        for (int z = std::max(0, z0); z <= std::min(nz - 1, z1); z++)
            for (int x = std::max(0, x0); x <= std::min(nx - 1, x1); x++) {
                real px = origin.x + x * spacing, pz = origin.z + z * spacing;
                real dx = (px - wx) / rx, dz = (pz - wz) / rz;
                real q = dx * dx + dz * dz;
                if (q >= 1) continue;
                real target = h[z * nx + x] - depth * (1 - q) * (1 - q);
                // stamp: keep the deepest impression at this sample
                real& cur = h[z * nx + x];
                real pressed = cur - depth * (1 - q) * (1 - q);
                if (pressed < cur) cur = pressed;
                (void)target;
            }
    }
    real sample(real wx, real wz) const {                  // bilinear height at world (x,z)
        real fx = (wx - origin.x) / spacing, fz = (wz - origin.z) / spacing;
        int x = (int)std::floor(fx), z = (int)std::floor(fz);
        if (x < 0) x = 0; if (x > nx - 2) x = nx - 2;
        if (z < 0) z = 0; if (z > nz - 2) z = nz - 2;
        real tx = fx - x, tz = fz - z;
        if (tx < 0) tx = 0; if (tx > 1) tx = 1; if (tz < 0) tz = 0; if (tz > 1) tz = 1;
        real h00 = h[z * nx + x], h10 = h[z * nx + x + 1], h01 = h[(z + 1) * nx + x], h11 = h[(z + 1) * nx + x + 1];
        return origin.y + (h00 * (1 - tx) + h10 * tx) * (1 - tz) + (h01 * (1 - tx) + h11 * tx) * tz;
    }
    Vector3 normal(real wx, real wz) const {               // via central differences
        real e = spacing * (real)0.5;
        Vector3 n(-(sample(wx + e, wz) - sample(wx - e, wz)) / (2 * e), 1,
                  -(sample(wx, wz + e) - sample(wx, wz - e)) / (2 * e));
        n.normalise(); return n;
    }
    unsigned sphereContact(const CollisionSphere& s, CollisionData* data) const {
        if (data->contactsLeft <= 0) return 0;
        Vector3 c = s.getAxis(3);
        real ground = sample(c.x, c.z);
        Vector3 n = normal(c.x, c.z);
        real pen = s.radius - (c.y - ground) * n.y;        // signed distance approx along normal
        if (pen <= 0) return 0;
        Contact* ct = data->contacts;
        ct->contactNormal = n; ct->penetration = pen;
        ct->contactPoint = c - n * s.radius;
        ct->setBodyData(s.body, nullptr, data->friction, data->restitution);
        data->addContacts(1); return 1;
    }
    unsigned capsuleContact(const CollisionCapsule& cap, CollisionData* data) const {
        Vector3 e[2]; cap.endpoints(e[0], e[1]); unsigned used = 0;
        for (int i = 0; i < 2 && data->hasMoreContacts(); i++) {
            real ground = sample(e[i].x, e[i].z); Vector3 n = normal(e[i].x, e[i].z);
            real pen = cap.radius - (e[i].y - ground) * n.y;
            if (pen <= 0) continue;
            Contact* ct = data->contacts;
            ct->contactNormal = n; ct->penetration = pen;
            ct->contactPoint = e[i] - n * cap.radius;
            ct->setBodyData(cap.body, nullptr, data->friction, data->restitution);
            data->addContacts(1); used++;
        }
        return used;
    }
    unsigned boxContact(const CollisionBox& box, CollisionData* data) const {
        static const real m[8][3] = {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
        unsigned used = 0;
        for (int i = 0; i < 8 && data->hasMoreContacts(); i++) {
            Vector3 v(m[i][0], m[i][1], m[i][2]); v.componentProductUpdate(box.halfSize);
            v = box.getTransform().transform(v);
            real ground = sample(v.x, v.z);
            if (v.y >= ground) continue;
            Vector3 n = normal(v.x, v.z);
            Contact* ct = data->contacts;
            ct->contactNormal = n; ct->penetration = ground - v.y;
            ct->contactPoint = v;
            ct->setBodyData(box.body, nullptr, data->friction, data->restitution);
            data->addContacts(1); used++;
        }
        return used;
    }
    bool raycast(const Vector3& origin_, const Vector3& dir, real maxDist, Vector3& hit) const {
        Vector3 d = dir.unit(); real step = spacing * (real)0.5;
        for (real t = 0; t < maxDist; t += step) {
            Vector3 p = origin_ + d * t;
            if (p.y <= sample(p.x, p.z)) {                 // bisect the last interval
                real lo = t - step > 0 ? t - step : 0, hi = t;
                for (int i = 0; i < 16; i++) { real mid = (lo + hi) / 2; Vector3 q = origin_ + d * mid;
                    if (q.y <= sample(q.x, q.z)) hi = mid; else lo = mid; }
                hit = origin_ + d * hi; return true;
            }
        }
        return false;
    }
};

} // namespace phys
