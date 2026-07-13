// Core linear-algebra types: Vector3, Quaternion, Matrix3 (3x3) and Matrix4
// (3x4 affine transform). Row-major storage. This is the mathematical
// foundation of the whole engine (Millington, Part I ch.2 and Part III ch.9).
#pragma once
#include "precision.h"

namespace phys {

class Vector3 {
public:
    real x, y, z;
private:
    real pad;                       // 4-word alignment, as in the book
public:
    Vector3() : x(0), y(0), z(0), pad(0) {}
    Vector3(real x, real y, real z) : x(x), y(y), z(z), pad(0) {}

    void clear() { x = y = z = 0; }
    void invert() { x = -x; y = -y; z = -z; }

    real magnitude() const { return real_sqrt(x * x + y * y + z * z); }
    real squareMagnitude() const { return x * x + y * y + z * z; }

    void normalise() { real l = magnitude(); if (l > 0) (*this) *= ((real)1) / l; }
    Vector3 unit() const { Vector3 r = *this; r.normalise(); return r; }

    void operator*=(real v) { x *= v; y *= v; z *= v; }
    Vector3 operator*(real v) const { return Vector3(x * v, y * v, z * v); }
    void operator+=(const Vector3& v) { x += v.x; y += v.y; z += v.z; }
    Vector3 operator+(const Vector3& v) const { return Vector3(x + v.x, y + v.y, z + v.z); }
    void operator-=(const Vector3& v) { x -= v.x; y -= v.y; z -= v.z; }
    Vector3 operator-(const Vector3& v) const { return Vector3(x - v.x, y - v.y, z - v.z); }
    Vector3 operator-() const { return Vector3(-x, -y, -z); }

    void addScaledVector(const Vector3& v, real t) { x += v.x * t; y += v.y * t; z += v.z * t; }

    // component (Hadamard) product
    Vector3 componentProduct(const Vector3& v) const { return Vector3(x * v.x, y * v.y, z * v.z); }
    void componentProductUpdate(const Vector3& v) { x *= v.x; y *= v.y; z *= v.z; }

    // scalar (dot) product
    real scalarProduct(const Vector3& v) const { return x * v.x + y * v.y + z * v.z; }
    real operator*(const Vector3& v) const { return x * v.x + y * v.y + z * v.z; }

    // vector (cross) product
    Vector3 vectorProduct(const Vector3& v) const {
        return Vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    Vector3 operator%(const Vector3& v) const { return vectorProduct(v); }
    void operator%=(const Vector3& v) { *this = vectorProduct(v); }

    real operator[](unsigned i) const { return i == 0 ? x : (i == 1 ? y : z); }
    real& operator[](unsigned i) { return i == 0 ? x : (i == 1 ? y : z); }

    bool operator==(const Vector3& o) const { return x == o.x && y == o.y && z == o.z; }

    // Build a right-handed orthonormal basis from a (normalised) and b.
    // On return a,b,c are mutually orthogonal unit vectors with a preserved.
    static void makeOrthonormalBasis(Vector3* a, Vector3* b, Vector3* c) {
        a->normalise();
        *c = (*a) % (*b);
        if (c->squareMagnitude() == 0) return;   // a and b were parallel
        c->normalise();
        *b = (*c) % (*a);
    }

    static const Vector3 GRAVITY;
    static const Vector3 UP;
    static const Vector3 X, Y, Z;
};

inline const Vector3 Vector3::GRAVITY(0, -9.81, 0);
inline const Vector3 Vector3::UP(0, 1, 0);
inline const Vector3 Vector3::X(1, 0, 0);
inline const Vector3 Vector3::Y(0, 1, 0);
inline const Vector3 Vector3::Z(0, 0, 1);

// -------------------------------------------------------------------- Quaternion
class Quaternion {
public:
    union {
        struct { real r, i, j, k; };
        real data[4];
    };
    Quaternion() : r(1), i(0), j(0), k(0) {}
    Quaternion(real r, real i, real j, real k) : r(r), i(i), j(j), k(k) {}

    void normalise() {
        real d = r * r + i * i + j * j + k * k;
        if (d < REAL_EPSILON) { r = 1; return; }
        d = ((real)1) / real_sqrt(d);
        r *= d; i *= d; j *= d; k *= d;
    }
    // Hamilton product (this = this * m)
    void operator*=(const Quaternion& m) {
        Quaternion q = *this;
        r = q.r * m.r - q.i * m.i - q.j * m.j - q.k * m.k;
        i = q.r * m.i + q.i * m.r + q.j * m.k - q.k * m.j;
        j = q.r * m.j + q.j * m.r + q.k * m.i - q.i * m.k;
        k = q.r * m.k + q.k * m.r + q.i * m.j - q.j * m.i;
    }
    // rotate this quaternion by a vector (scaled orientation update helper)
    void rotateByVector(const Vector3& v) { Quaternion q(0, v.x, v.y, v.z); (*this) *= q; }
    // q += (1/2) * scale * omega * q     — the quaternion derivative from angular velocity
    void addScaledVector(const Vector3& v, real scale) {
        Quaternion q(0, v.x * scale, v.y * scale, v.z * scale);
        q *= *this;
        r += q.r * ((real)0.5); i += q.i * ((real)0.5);
        j += q.j * ((real)0.5); k += q.k * ((real)0.5);
    }
};

// ---------------------------------------------------------------------- Matrix3
class Matrix3 {
public:
    real data[9];                    // row-major 3x3
    Matrix3() { for (int n = 0; n < 9; n++) data[n] = 0; }
    Matrix3(real c0, real c1, real c2, real c3, real c4, real c5, real c6, real c7, real c8) {
        data[0] = c0; data[1] = c1; data[2] = c2; data[3] = c3; data[4] = c4;
        data[5] = c5; data[6] = c6; data[7] = c7; data[8] = c8;
    }
    void setDiagonal(real a, real b, real c) { setInertiaTensorCoeffs(a, b, c); }
    void setInertiaTensorCoeffs(real ix, real iy, real iz, real ixy = 0, real ixz = 0, real iyz = 0) {
        data[0] = ix;  data[1] = -ixy; data[2] = -ixz;
        data[3] = -ixy; data[4] = iy;  data[5] = -iyz;
        data[6] = -ixz; data[7] = -iyz; data[8] = iz;
    }
    // inertia tensor of a rectangular block of given half-sizes and mass
    void setBlockInertiaTensor(const Vector3& halfSizes, real mass) {
        Vector3 sq = halfSizes.componentProduct(halfSizes);
        setInertiaTensorCoeffs(((real)0.3) * mass * (sq.y + sq.z),
                               ((real)0.3) * mass * (sq.x + sq.z),
                               ((real)0.3) * mass * (sq.x + sq.y));
    }
    void setSkewSymmetric(const Vector3& v) {
        data[0] = data[4] = data[8] = 0;
        data[1] = -v.z; data[2] = v.y;
        data[3] = v.z;  data[5] = -v.x;
        data[6] = -v.y; data[7] = v.x;
    }
    void setComponents(const Vector3& a, const Vector3& b, const Vector3& c) {
        data[0] = a.x; data[1] = b.x; data[2] = c.x;
        data[3] = a.y; data[4] = b.y; data[5] = c.y;
        data[6] = a.z; data[7] = b.z; data[8] = c.z;
    }
    Vector3 operator*(const Vector3& v) const {
        return Vector3(v.x * data[0] + v.y * data[1] + v.z * data[2],
                       v.x * data[3] + v.y * data[4] + v.z * data[5],
                       v.x * data[6] + v.y * data[7] + v.z * data[8]);
    }
    Vector3 transform(const Vector3& v) const { return (*this) * v; }
    Vector3 transformTranspose(const Vector3& v) const {
        return Vector3(v.x * data[0] + v.y * data[3] + v.z * data[6],
                       v.x * data[1] + v.y * data[4] + v.z * data[7],
                       v.x * data[2] + v.y * data[5] + v.z * data[8]);
    }
    Matrix3 operator*(const Matrix3& o) const {
        return Matrix3(
            data[0]*o.data[0]+data[1]*o.data[3]+data[2]*o.data[6],
            data[0]*o.data[1]+data[1]*o.data[4]+data[2]*o.data[7],
            data[0]*o.data[2]+data[1]*o.data[5]+data[2]*o.data[8],
            data[3]*o.data[0]+data[4]*o.data[3]+data[5]*o.data[6],
            data[3]*o.data[1]+data[4]*o.data[4]+data[5]*o.data[7],
            data[3]*o.data[2]+data[4]*o.data[5]+data[5]*o.data[8],
            data[6]*o.data[0]+data[7]*o.data[3]+data[8]*o.data[6],
            data[6]*o.data[1]+data[7]*o.data[4]+data[8]*o.data[7],
            data[6]*o.data[2]+data[7]*o.data[5]+data[8]*o.data[8]);
    }
    void operator*=(const Matrix3& o) { *this = (*this) * o; }
    void operator*=(real s) { for (int n = 0; n < 9; n++) data[n] *= s; }
    void operator+=(const Matrix3& o) { for (int n = 0; n < 9; n++) data[n] += o.data[n]; }

    Vector3 getRowVector(int i) const { return Vector3(data[i * 3], data[i * 3 + 1], data[i * 3 + 2]); }
    Vector3 getAxisVector(int i) const { return Vector3(data[i], data[i + 3], data[i + 6]); }

    void setInverse(const Matrix3& m) {
        real a = m.data[0], b = m.data[1], c = m.data[2];
        real d = m.data[3], e = m.data[4], f = m.data[5];
        real g = m.data[6], h = m.data[7], i = m.data[8];
        real det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
        if (det == 0) return;
        real id = ((real)1) / det;
        data[0] = (e * i - f * h) * id;  data[1] = (c * h - b * i) * id;  data[2] = (b * f - c * e) * id;
        data[3] = (f * g - d * i) * id;  data[4] = (a * i - c * g) * id;  data[5] = (c * d - a * f) * id;
        data[6] = (d * h - e * g) * id;  data[7] = (b * g - a * h) * id;  data[8] = (a * e - b * d) * id;
    }
    Matrix3 inverse() const { Matrix3 r; r.setInverse(*this); return r; }
    void invert() { setInverse(Matrix3(*this)); }

    void setTranspose(const Matrix3& m) {
        data[0] = m.data[0]; data[1] = m.data[3]; data[2] = m.data[6];
        data[3] = m.data[1]; data[4] = m.data[4]; data[5] = m.data[7];
        data[6] = m.data[2]; data[7] = m.data[5]; data[8] = m.data[8];
    }
    Matrix3 transpose() const { Matrix3 r; r.setTranspose(*this); return r; }

    void setOrientation(const Quaternion& q) {
        data[0] = 1 - 2 * (q.j * q.j + q.k * q.k);
        data[1] = 2 * (q.i * q.j - q.k * q.r);
        data[2] = 2 * (q.i * q.k + q.j * q.r);
        data[3] = 2 * (q.i * q.j + q.k * q.r);
        data[4] = 1 - 2 * (q.i * q.i + q.k * q.k);
        data[5] = 2 * (q.j * q.k - q.i * q.r);
        data[6] = 2 * (q.i * q.k - q.j * q.r);
        data[7] = 2 * (q.j * q.k + q.i * q.r);
        data[8] = 1 - 2 * (q.i * q.i + q.j * q.j);
    }
    static Matrix3 linearInterpolate(const Matrix3& a, const Matrix3& b, real prop) {
        Matrix3 r; for (int n = 0; n < 9; n++) r.data[n] = a.data[n] * (1 - prop) + b.data[n] * prop; return r;
    }
};

// ---------------------------------------------------------------------- Matrix4
// 3x4 affine transform (rotation/scale in the 3x3 block, translation in col 3).
class Matrix4 {
public:
    real data[12];
    Matrix4() { for (int n = 0; n < 12; n++) data[n] = 0; data[0] = data[5] = data[10] = 1; }

    Vector3 operator*(const Vector3& v) const {           // transform point
        return Vector3(v.x * data[0] + v.y * data[1] + v.z * data[2] + data[3],
                       v.x * data[4] + v.y * data[5] + v.z * data[6] + data[7],
                       v.x * data[8] + v.y * data[9] + v.z * data[10] + data[11]);
    }
    Vector3 transform(const Vector3& v) const { return (*this) * v; }
    Vector3 transformDirection(const Vector3& v) const {
        return Vector3(v.x * data[0] + v.y * data[1] + v.z * data[2],
                       v.x * data[4] + v.y * data[5] + v.z * data[6],
                       v.x * data[8] + v.y * data[9] + v.z * data[10]);
    }
    // inverse-transform assuming an orthonormal 3x3 block (rigid transform)
    Vector3 transformInverse(const Vector3& v) const {
        Vector3 t = v; t.x -= data[3]; t.y -= data[7]; t.z -= data[11];
        return Vector3(t.x * data[0] + t.y * data[4] + t.z * data[8],
                       t.x * data[1] + t.y * data[5] + t.z * data[9],
                       t.x * data[2] + t.y * data[6] + t.z * data[10]);
    }
    Vector3 transformInverseDirection(const Vector3& v) const {
        return Vector3(v.x * data[0] + v.y * data[4] + v.z * data[8],
                       v.x * data[1] + v.y * data[5] + v.z * data[9],
                       v.x * data[2] + v.y * data[6] + v.z * data[10]);
    }
    Vector3 getAxisVector(int i) const { return Vector3(data[i], data[i + 4], data[i + 8]); }

    Matrix4 operator*(const Matrix4& o) const {           // compose affine transforms
        Matrix4 r;
        r.data[0]  = o.data[0]*data[0]  + o.data[4]*data[1]  + o.data[8]*data[2];
        r.data[4]  = o.data[0]*data[4]  + o.data[4]*data[5]  + o.data[8]*data[6];
        r.data[8]  = o.data[0]*data[8]  + o.data[4]*data[9]  + o.data[8]*data[10];
        r.data[1]  = o.data[1]*data[0]  + o.data[5]*data[1]  + o.data[9]*data[2];
        r.data[5]  = o.data[1]*data[4]  + o.data[5]*data[5]  + o.data[9]*data[6];
        r.data[9]  = o.data[1]*data[8]  + o.data[5]*data[9]  + o.data[9]*data[10];
        r.data[2]  = o.data[2]*data[0]  + o.data[6]*data[1]  + o.data[10]*data[2];
        r.data[6]  = o.data[2]*data[4]  + o.data[6]*data[5]  + o.data[10]*data[6];
        r.data[10] = o.data[2]*data[8]  + o.data[6]*data[9]  + o.data[10]*data[10];
        r.data[3]  = o.data[3]*data[0]  + o.data[7]*data[1]  + o.data[11]*data[2]  + data[3];
        r.data[7]  = o.data[3]*data[4]  + o.data[7]*data[5]  + o.data[11]*data[6]  + data[7];
        r.data[11] = o.data[3]*data[8]  + o.data[7]*data[9]  + o.data[11]*data[10] + data[11];
        return r;
    }

    real getDeterminant() const {
        return -data[8] * data[5] * data[2] + data[4] * data[9] * data[2] + data[8] * data[1] * data[6]
               - data[0] * data[9] * data[6] - data[4] * data[1] * data[10] + data[0] * data[5] * data[10];
    }
    void setInverse(const Matrix4& m) {
        real det = m.getDeterminant();
        if (det == 0) return;
        det = ((real)1) / det;
        data[0]  = (-m.data[9]*m.data[6]  + m.data[5]*m.data[10]) * det;
        data[4]  = ( m.data[8]*m.data[6]  - m.data[4]*m.data[10]) * det;
        data[8]  = (-m.data[8]*m.data[5]  + m.data[4]*m.data[9] ) * det;
        data[1]  = ( m.data[9]*m.data[2]  - m.data[1]*m.data[10]) * det;
        data[5]  = (-m.data[8]*m.data[2]  + m.data[0]*m.data[10]) * det;
        data[9]  = ( m.data[8]*m.data[1]  - m.data[0]*m.data[9] ) * det;
        data[2]  = (-m.data[5]*m.data[2]  + m.data[1]*m.data[6] ) * det;
        data[6]  = ( m.data[4]*m.data[2]  - m.data[0]*m.data[6] ) * det;
        data[10] = (-m.data[4]*m.data[1]  + m.data[0]*m.data[5] ) * det;
        data[3]  = ( m.data[9]*m.data[6]*m.data[3]  - m.data[5]*m.data[10]*m.data[3]
                   - m.data[9]*m.data[2]*m.data[7]  + m.data[1]*m.data[10]*m.data[7]
                   + m.data[5]*m.data[2]*m.data[11] - m.data[1]*m.data[6]*m.data[11]) * det;
        data[7]  = (-m.data[8]*m.data[6]*m.data[3]  + m.data[4]*m.data[10]*m.data[3]
                   + m.data[8]*m.data[2]*m.data[7]  - m.data[0]*m.data[10]*m.data[7]
                   - m.data[4]*m.data[2]*m.data[11] + m.data[0]*m.data[6]*m.data[11]) * det;
        data[11] = ( m.data[8]*m.data[5]*m.data[3]  - m.data[4]*m.data[9]*m.data[3]
                   - m.data[8]*m.data[1]*m.data[7]  + m.data[0]*m.data[9]*m.data[7]
                   + m.data[4]*m.data[1]*m.data[11] - m.data[0]*m.data[5]*m.data[11]) * det;
    }
    Matrix4 inverse() const { Matrix4 r; r.setInverse(*this); return r; }
    void invert() { setInverse(Matrix4(*this)); }

    void setOrientationAndPos(const Quaternion& q, const Vector3& pos) {
        data[0] = 1 - 2 * (q.j * q.j + q.k * q.k);
        data[1] = 2 * (q.i * q.j - q.k * q.r);
        data[2] = 2 * (q.i * q.k + q.j * q.r);
        data[3] = pos.x;
        data[4] = 2 * (q.i * q.j + q.k * q.r);
        data[5] = 1 - 2 * (q.i * q.i + q.k * q.k);
        data[6] = 2 * (q.j * q.k - q.i * q.r);
        data[7] = pos.y;
        data[8] = 2 * (q.i * q.k - q.j * q.r);
        data[9] = 2 * (q.j * q.k + q.i * q.r);
        data[10] = 1 - 2 * (q.i * q.i + q.j * q.j);
        data[11] = pos.z;
    }
};

} // namespace phys
