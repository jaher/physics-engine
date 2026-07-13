// Volumetric deformables + PBF fluids (phys/softbody.h):
//   • polar decomposition returns a rotation for a rotation and identity for a stretch
//   • a soft box dropped on the ground squashes then rebounds, conserving volume
//   • the self-collision helper separates overlapping particles to ≥ 2·radius
//   • a PBF pool holds its rest density after settling
#include "phys/softbody.h"
#include "check.h"
#include <cmath>
#include <vector>
using namespace phys;

// Rotation matrix about a (unit) axis by angle, via the engine's quaternion path.
static Matrix3 rotationMatrix(Vector3 axis, real ang) {
    axis.normalise();
    real s = real_sin(ang / 2), c = real_cos(ang / 2);
    Quaternion q(c, axis.x * s, axis.y * s, axis.z * s); q.normalise();
    Matrix3 m; m.setOrientation(q); return m;
}
static bool matNear(const Matrix3& a, const Matrix3& b, real tol) {
    for (int i = 0; i < 9; ++i) if (real_abs(a.data[i] - b.data[i]) > tol) return false;
    return true;
}

int main() {
    // A) polar decomposition -------------------------------------------------
    {
        Matrix3 R = rotationMatrix(Vector3(1, 2, -1), (real)0.7);
        Matrix3 Rp, Sp; polarDecompose(R, Rp, Sp);
        CHECK(matNear(Rp, R, 1e-6));                 // rotation → itself
        CHECK(matNear(Sp, Matrix3(1,0,0, 0,1,0, 0,0,1), 1e-6));  // no stretch (S ≈ I)

        Matrix3 stretch(1.4,0,0, 0,0.7,0, 0,0,1.2);  // symmetric → pure stretch
        Matrix3 Rs, Ss; polarDecompose(stretch, Rs, Ss);
        CHECK(matNear(Rs, Matrix3(1,0,0, 0,1,0, 0,0,1), 1e-6));  // no rotation
        CHECK(matNear(Ss, stretch, 1e-6));           // stretch recovered

        // general F = R·S → polar recovers R
        Matrix3 F = R * stretch;
        CHECK(matNear(polarRotation(F), R, 1e-5));
    }

    // B) soft box drop: squash + rebound + volume conservation ---------------
    {
        // 0.24 m cube of a soft solid, its base 0.12 m above the ground.
        SoftBody sb = SoftBody::makeBox(3, 3, 3, (real)0.08, /*E*/15000, /*nu*/0.42,
                                        /*rho*/600, Vector3(-0.12, 0.12, -0.12));
        sb.groundY = 0; sb.restitution = (real)0.4; sb.friction = (real)0.25;
        real rest = sb.restVolume();
        CHECK(rest > 0);
        CHECK((int)sb.tets.size() == 3 * 3 * 3 * 6);   // 6 tets per cell

        real dt = 1e-3;
        real restH = 3 * (real)0.08;                    // rest vertical extent
        real minH = restH, minCom = 1e9;
        real volLo = rest, volHi = rest;
        bool finite = true, aboveGround = true, rebounded = false;
        for (int s = 0; s < 2500; ++s) {
            sb.step(dt);
            // geometry
            real lo = 1e9, hi = -1e9, com = 0;
            for (auto& n : sb.nodes) { lo = std::min(lo, n.x.y); hi = std::max(hi, n.x.y); com += n.x.y;
                if (!std::isfinite(n.x.y)) finite = false; }
            com /= sb.nodes.size();
            real H = hi - lo;
            minH = std::min(minH, H);
            if (com < minCom) minCom = com;
            else if (minCom < 1e8 && com > minCom + (real)0.001) rebounded = true;  // rose after its low point
            if (lo < -(real)1e-2) aboveGround = false;
            real v = sb.totalVolume();
            volLo = std::min(volLo, v); volHi = std::max(volHi, v);
        }
        CHECK(finite);
        CHECK(aboveGround);                              // never tunnels through the floor
        CHECK(minH < restH * (real)0.99);               // visibly squashed
        CHECK(rebounded);                                // COM rose again after impact
        CHECK(volLo > rest * (real)0.80);               // volume within ~20% of rest…
        CHECK(volHi < rest * (real)1.20);               // …throughout the whole drop
    }

    // C) self-collision: two overlapping particles pushed to ≥ 2·radius -------
    {
        real radius = (real)0.1;
        std::vector<Vector3> pos = { Vector3(0, 0, 0), Vector3(0.05, 0, 0) };  // 0.05 < 0.2 apart
        resolveSelfCollisions(pos, radius);
        real sep = (pos[0] - pos[1]).magnitude();
        CHECK(sep >= 2 * radius - (real)1e-9);

        // a tight clump of many particles ends with no remaining deep overlaps
        std::vector<Vector3> clump;
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j)
            clump.push_back(Vector3(i * 0.03, j * 0.03, 0));   // spacing 0.03 << 0.2
        resolveSelfCollisions(clump, radius, 8);
        real minSep = 1e9;
        for (size_t a = 0; a < clump.size(); ++a)
            for (size_t b = a + 1; b < clump.size(); ++b)
                minSep = std::min(minSep, (double)(clump[a] - clump[b]).magnitude());
        CHECK(minSep > 2 * radius * 0.85);              // no pair left deeply overlapped
    }

    // D) PBF keeps mean density within ~10% of rest after settling -----------
    {
        PBFFluid f; f.h = (real)0.1; f.solverIters = 5; f.xsph = (real)0.02;
        Vector3 cmin(-0.2, 0, -0.2), cmax(0.2, 0.4, 0.2);
        f.build(cmin, cmax, Vector3(-0.18, 0.03, -0.18), Vector3(0.18, 0.38, 0.18), (real)0.05);
        int n0 = f.nParticles();
        CHECK(n0 > 200);
        CHECK(f.rho0 > 0);
        real rest = f.rho0;
        for (int s = 0; s < 300; ++s) f.step(2e-3);
        CHECK(f.nParticles() == n0);                    // no particles created/destroyed
        bool finite = true, inside = true;
        for (int i = 0; i < n0; ++i) {
            if (!std::isfinite(f.pos[i].squareMagnitude())) finite = false;
            for (int a = 0; a < 3; ++a)
                if (f.pos[i][a] < cmin[a] - 1e-3 || f.pos[i][a] > cmax[a] + 1e-3) inside = false;
        }
        CHECK(finite);
        CHECK(inside);
        real mean = f.meanDensity();
        CHECK(mean > rest * (real)0.90);
        CHECK(mean < rest * (real)1.10);
    }

    return test::report("softbody");
}
