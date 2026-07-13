#include "phys/core.h"
#include "check.h"
using namespace phys;

static bool vnear(const Vector3& a, const Vector3& b, double t = 1e-6) {
    return test::near(a.x, b.x, t) && test::near(a.y, b.y, t) && test::near(a.z, b.z, t);
}

int main() {
    // dot / cross
    Vector3 a(1, 2, 3), b(4, 5, 6);
    CHECK_NEAR(a * b, 32.0, 1e-9);
    CHECK(vnear(a % b, Vector3(-3, 6, -3)));
    CHECK(vnear(Vector3::X % Vector3::Y, Vector3::Z));      // right-handed

    // magnitude / normalise
    Vector3 v(3, 4, 0);
    CHECK_NEAR(v.magnitude(), 5.0, 1e-9);
    CHECK_NEAR(v.unit().magnitude(), 1.0, 1e-9);

    // Matrix3 inverse: M * M^-1 = I
    Matrix3 M(1, 2, 3, 0, 1, 4, 5, 6, 0);
    Matrix3 I = M * M.inverse();
    CHECK_NEAR(I.data[0], 1, 1e-9); CHECK_NEAR(I.data[4], 1, 1e-9); CHECK_NEAR(I.data[8], 1, 1e-9);
    CHECK_NEAR(I.data[1], 0, 1e-9); CHECK_NEAR(I.data[5], 0, 1e-9); CHECK_NEAR(I.data[6], 0, 1e-9);

    // Quaternion → matrix: 90° about Z maps X→Y
    Quaternion q(std::cos(real_pi / 4), 0, 0, std::sin(real_pi / 4));
    q.normalise();
    Matrix3 R; R.setOrientation(q);
    CHECK(vnear(R * Vector3::X, Vector3::Y, 1e-9));
    CHECK(vnear(R.transpose() * (R * a), a, 1e-9));         // R orthonormal

    // Matrix4 transform + inverse round trip (rotate 90° Z, translate)
    Matrix4 T; T.setOrientationAndPos(q, Vector3(10, 0, 0));
    Vector3 p(1, 0, 0), wp = T * p;
    CHECK(vnear(wp, Vector3(10, 1, 0), 1e-9));
    CHECK(vnear(T.transformInverse(wp), p, 1e-9));
    CHECK(vnear(T.inverse() * wp, p, 1e-9));

    // quaternion integration: many small Euler steps of spin about Y converge to
    // a true 90° rotation, mapping Z→X (validates addScaledVector direction).
    Quaternion o; Vector3 omega(0, 1, 0);
    real dt = (real_pi / 2) / 4000;
    for (int s = 0; s < 4000; s++) { o.addScaledVector(omega, dt); o.normalise(); }
    Matrix3 RO; RO.setOrientation(o);
    CHECK(vnear(RO * Vector3::Z, Vector3::X, 1e-3));

    return test::report("core");
}
