// Reduced-coordinate articulated-body dynamics — Featherstone's Articulated-Body
// Algorithm (ABA), the O(n) robotics standard the maximal-coordinate joint
// solvers (joints.h/joints2.h) lack. A serial/tree chain of 1-DOF joints
// (REVOLUTE or PRISMATIC) is simulated in *joint* coordinates: applied joint
// torques + gravity → joint accelerations, then semi-implicit-Euler integrated.
//
// Spatial (6D) formulation, Featherstone "Rigid Body Dynamics Algorithms" (2008)
// ch.5-7. All spatial vectors are expressed in the WORLD inertial frame and
// taken about the world origin, so every parent→child Plücker transform is the
// identity — the three ABA sweeps then reduce to plain adds. This is exact: the
// only frame-dependent term, the velocity product c = v ×ₘ (S q̇), is unchanged
// because (S q̇) ×ₘ (S q̇) = 0, so parent- and body-frame forms coincide.
//
// Gravity is applied by Featherstone's base-acceleration trick: the (fixed) base
// is given spatial acceleration a₀ = [0; -g], making every body feel its weight
// without an explicit external force. A single released link then reproduces the
// compound-pendulum α = m g l / I_pivot; see tests/test_articulation.cpp.
#pragma once
#include "core.h"
#include <vector>

namespace phys {

// ------------------------------------------------------------ spatial 6-vector
// Motion vectors store [ang = ω; lin = v_O]; force vectors store [ang = n; lin = f].
struct SpatialVector {
    Vector3 ang, lin;
    SpatialVector() {}
    SpatialVector(const Vector3& a, const Vector3& l) : ang(a), lin(l) {}
    SpatialVector operator+(const SpatialVector& o) const { return SpatialVector(ang + o.ang, lin + o.lin); }
    SpatialVector operator-(const SpatialVector& o) const { return SpatialVector(ang - o.ang, lin - o.lin); }
    SpatialVector operator*(real s) const { return SpatialVector(ang * s, lin * s); }
    void operator+=(const SpatialVector& o) { ang += o.ang; lin += o.lin; }
};

// ------------------------------------------------------------ spatial 6x6 inertia
// Stored as four 3x3 blocks; force = M·motion.  Rigid-body inertias start of the
// rigid form and accumulate into general symmetric articulated inertias in ABA.
struct SpatialMatrix {
    Matrix3 m11, m12, m21, m22;
    SpatialVector operator*(const SpatialVector& v) const {
        return SpatialVector(m11 * v.ang + m12 * v.lin, m21 * v.ang + m22 * v.lin);
    }
    void operator+=(const SpatialMatrix& o) { m11 += o.m11; m12 += o.m12; m21 += o.m21; m22 += o.m22; }
};

enum class JointType { REVOLUTE, PRISMATIC };

class Articulation {
public:
    // Add a 1-DOF link hanging off `parent` (-1 = fixed ground). `axis` is the
    // joint axis and `jointPos` the joint anchor, both in the parent's frame at
    // q=0; `comOffset`/`inertia` (about the COM) are in the link's own frame.
    int addLink(int parent, JointType type, Vector3 axis, real mass,
                Matrix3 inertia, Vector3 comOffset, Vector3 jointPos) {
        Link L;
        L.parent = parent; L.type = type;
        L.axis = axis.unit(); L.mass = mass; L.inertiaBody = inertia;
        L.comOffset = comOffset; L.jointPos = jointPos;
        L.q = L.qd = L.tau = L.qdd = 0;
        links.push_back(L);
        return (int)links.size() - 1;
    }

    void setJointTorque(int link, real tau) { links[link].tau = tau; }
    void setGravity(Vector3 g) { gravity = g; }
    void setJointState(int link, real q, real qd) { links[link].q = q; links[link].qd = qd; }

    real q(int link)     const { return links[link].q; }
    real qd(int link)    const { return links[link].qd; }
    real qddot(int link) const { return links[link].qdd; }   // last computed acceleration
    int  numLinks()      const { return (int)links.size(); }

    // World transform of a link's frame (origin at its joint anchor, orientation R).
    Matrix4 linkWorld(int link) {
        computeKinematics();
        const Link& L = links[link];
        Matrix4 m;
        m.data[0] = L.R.data[0]; m.data[1] = L.R.data[1]; m.data[2]  = L.R.data[2]; m.data[3]  = L.A.x;
        m.data[4] = L.R.data[3]; m.data[5] = L.R.data[4]; m.data[6]  = L.R.data[5]; m.data[7]  = L.A.y;
        m.data[8] = L.R.data[6]; m.data[9] = L.R.data[7]; m.data[10] = L.R.data[8]; m.data[11] = L.A.z;
        return m;
    }
    Vector3 comWorld(int link)    { computeKinematics(); return links[link].C; }
    Vector3 linkOmega(int link)   { computeKinematics(); computeVelocities(); return links[link].v.ang; }

    // Total mechanical energy (KE + PE) — the invariant the tests watch drift on.
    real totalEnergy() {
        computeKinematics();
        computeVelocities();
        real E = 0;
        for (const Link& L : links) {
            Vector3 w = L.v.ang;
            Vector3 vG = L.v.lin + w % L.C;                    // COM velocity = v_O + ω×c
            real ke = ((real)0.5) * L.mass * vG.scalarProduct(vG)
                    + ((real)0.5) * w.scalarProduct(L.Iworld * w);
            real pe = -L.mass * gravity.scalarProduct(L.C);    // U = -m g·x
            E += ke + pe;
        }
        return E;
    }

    // One ABA forward-dynamics solve + semi-implicit-Euler integration step.
    void step(real dt) {
        int n = (int)links.size();
        computeKinematics();
        computeVelocities();

        // Pass 1 (base→tip): velocity-product bias c, articulated inertia/bias seed.
        for (int i = 0; i < n; i++) {
            Link& L = links[i];
            SpatialVector vJ = jointMotion(L, L.qd);
            L.c  = crm(L.v, vJ);
            L.IA = L.Irb;
            L.pA = crf(L.v, L.Irb * L.v);                       // v ×f (I v) — gyroscopic bias
        }

        // Pass 2 (tip→base): articulated inertias/forces fold onto parents.
        for (int i = n - 1; i >= 0; i--) {
            Link& L = links[i];
            L.U = L.IA * L.S;
            L.D = dot(L.U, L.S);                                // S^T I^A S (scalar)
            L.u = L.tau - dot(L.pA, L.S);
            if (L.parent >= 0) {
                real invD = ((real)1) / L.D;
                SpatialMatrix Ia = subRank1(L.IA, L.U, invD);   // I^A − U D⁻¹ Uᵀ
                SpatialVector pa = L.pA + Ia * L.c + L.U * (L.u * invD);
                links[L.parent].IA += Ia;
                links[L.parent].pA += pa;
            }
        }

        // Pass 3 (base→tip): resolve accelerations. Fixed base carries a₀ = [0;-g].
        SpatialVector a0(Vector3(0, 0, 0), -gravity);
        for (int i = 0; i < n; i++) {
            Link& L = links[i];
            SpatialVector aPar = (L.parent >= 0) ? links[L.parent].a : a0;
            SpatialVector aPrime = aPar + L.c;
            L.qdd = (L.u - dot(L.U, aPrime)) / L.D;
            L.a = aPrime + jointMotion(L, L.qdd);
        }

        // Semi-implicit (symplectic) Euler: velocity first, then position.
        for (int i = 0; i < n; i++) {
            Link& L = links[i];
            L.qd += L.qdd * dt;
            L.q  += L.qd  * dt;
        }
    }

private:
    struct Link {
        int parent; JointType type;
        Vector3 axis, comOffset, jointPos; real mass; Matrix3 inertiaBody;
        real q, qd, tau, qdd;
        // cached world kinematics
        Matrix3 R, Iworld; Vector3 A, C, axisW;
        SpatialVector S;                 // joint motion subspace (world, about origin)
        SpatialMatrix Irb;               // rigid-body spatial inertia
        // ABA scratch
        SpatialVector v, c, pA, U, a; SpatialMatrix IA; real D, u;
    };
    std::vector<Link> links;
    Vector3 gravity = Vector3(0, (real)-9.81, 0);

    // --- spatial algebra -----------------------------------------------------
    static Vector3 mul(const Matrix3& m, const Vector3& v) { return m * v; }
    // motion×motion → motion
    static SpatialVector crm(const SpatialVector& v, const SpatialVector& m) {
        return SpatialVector(v.ang % m.ang, (v.ang % m.lin) + (v.lin % m.ang));
    }
    // motion×force → force
    static SpatialVector crf(const SpatialVector& v, const SpatialVector& f) {
        return SpatialVector((v.ang % f.ang) + (v.lin % f.lin), v.ang % f.lin);
    }
    static real dot(const SpatialVector& a, const SpatialVector& b) {
        return a.ang.scalarProduct(b.ang) + a.lin.scalarProduct(b.lin);
    }
    static Matrix3 outer(const Vector3& a, const Vector3& b) {
        return Matrix3(a.x*b.x, a.x*b.y, a.x*b.z,
                       a.y*b.x, a.y*b.y, a.y*b.z,
                       a.z*b.x, a.z*b.y, a.z*b.z);
    }
    // I^A − (1/D) U Uᵀ  (rank-1 downdate, block by block)
    static SpatialMatrix subRank1(const SpatialMatrix& I, const SpatialVector& U, real invD) {
        SpatialMatrix r = I;
        subScaled(r.m11, outer(U.ang, U.ang), invD);
        subScaled(r.m12, outer(U.ang, U.lin), invD);
        subScaled(r.m21, outer(U.lin, U.ang), invD);
        subScaled(r.m22, outer(U.lin, U.lin), invD);
        return r;
    }
    static void subScaled(Matrix3& a, const Matrix3& b, real s) {
        for (int k = 0; k < 9; k++) a.data[k] -= b.data[k] * s;
    }

    // Spatial motion produced by unit joint rate, scaled by `rate`.
    static SpatialVector jointMotion(const Link& L, real rate) {
        return SpatialVector(L.S.ang * rate, L.S.lin * rate);
    }

    // Axis-angle → rotation matrix via the core Quaternion.
    static Matrix3 rotAxis(const Vector3& axis, real angle) {
        real h = angle * (real)0.5, s = real_sin(h);
        Quaternion qq(real_cos(h), axis.x * s, axis.y * s, axis.z * s);
        Matrix3 R; R.setOrientation(qq); return R;
    }

    // Forward kinematics + per-link spatial inertia / joint subspace (world frame).
    void computeKinematics() {
        static const Matrix3 I3(1,0,0, 0,1,0, 0,0,1);
        for (Link& L : links) {
            Matrix3 Rp = (L.parent >= 0) ? links[L.parent].R : I3;
            Vector3 Ap = (L.parent >= 0) ? links[L.parent].A : Vector3();
            L.axisW = Rp * L.axis;                              // axis in world
            Vector3 anchor = Ap + Rp * L.jointPos;
            if (L.type == JointType::REVOLUTE) {
                L.R = Rp * rotAxis(L.axis, L.q);
                L.A = anchor;
                L.S = SpatialVector(L.axisW, L.A % L.axisW);    // [a; p×a]
            } else { // PRISMATIC
                L.R = Rp;
                L.A = anchor + L.axisW * L.q;
                L.S = SpatialVector(Vector3(), L.axisW);        // [0; a]
            }
            L.C = L.A + L.R * L.comOffset;                      // COM world position
            L.Iworld = L.R * L.inertiaBody * L.R.transpose();   // COM inertia in world

            // Rigid-body spatial inertia about the world origin (Featherstone 2.63):
            //   [ Ī − m c̃²   m c̃ ; −m c̃   m 1 ]
            Matrix3 cx; cx.setSkewSymmetric(L.C);
            Matrix3 cx2 = cx * cx;
            Matrix3 m11 = L.Iworld;
            for (int k = 0; k < 9; k++) m11.data[k] -= L.mass * cx2.data[k];
            Matrix3 mcx = cx; mcx *= L.mass;
            L.Irb.m11 = m11;
            L.Irb.m12 = mcx;
            L.Irb.m21 = mcx; L.Irb.m21 *= (real)-1;
            L.Irb.m22 = Matrix3(L.mass,0,0, 0,L.mass,0, 0,0,L.mass);
        }
    }

    // Spatial velocities: v_i = v_parent + S_i q̇_i  (identity transform).
    void computeVelocities() {
        for (Link& L : links) {
            SpatialVector vp = (L.parent >= 0) ? links[L.parent].v : SpatialVector(Vector3(), Vector3());
            L.v = vp + jointMotion(L, L.qd);
        }
    }
};

} // namespace phys
