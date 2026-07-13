#include "phys/fgen.h"
#include "check.h"
using namespace phys;

static RigidBody makeBox(real mass, Vector3 half) {
    RigidBody b; b.setMass(mass);
    Matrix3 I; I.setBlockInertiaTensor(half, mass);
    b.setInertiaTensor(I);
    b.setDamping(1, 1); b.setCanSleep(false);
    b.calculateDerivedData();
    return b;
}

int main() {
    // A) Free fall: same as a particle.
    {
        RigidBody b = makeBox(1, Vector3(0.5, 0.5, 0.5));
        b.setAcceleration(0, -9.81, 0); b.calculateDerivedData();
        for (int i = 0; i < 1000; i++) b.integrate(0.001);
        CHECK_NEAR(b.getVelocity().y, -9.81, 1e-6);
        CHECK_NEAR(b.getPosition().y, -0.5 * 9.81, 2e-2);
    }

    // B) Torque produces angular acceleration Iinv·τ.
    {
        RigidBody b; b.setMass(1);
        Matrix3 I; I.setInertiaTensorCoeffs(2, 3, 4); b.setInertiaTensor(I);
        b.setDamping(1, 1); b.setCanSleep(false); b.calculateDerivedData();
        b.addTorque(Vector3(1, 0, 0));
        b.integrate(0.01);
        CHECK_NEAR(b.getRotation().x, 0.5 * 0.01, 1e-6);          // Iinv_xx = 1/2
        CHECK_NEAR(b.getRotation().y, 0.0, 1e-9);
        CHECK_NEAR(b.getRotation().z, 0.0, 1e-9);
    }

    // C) A body spinning about a principal axis keeps that spin (torque-free).
    {
        RigidBody b; b.setMass(1);
        Matrix3 I; I.setInertiaTensorCoeffs(2, 3, 4); b.setInertiaTensor(I);
        b.setDamping(1, 1); b.setCanSleep(false); b.setRotation(0, 5, 0); b.calculateDerivedData();
        for (int i = 0; i < 2000; i++) b.integrate(0.001);
        CHECK_NEAR(b.getRotation().x, 0.0, 1e-3);
        CHECK_NEAR(b.getRotation().y, 5.0, 1e-3);
        CHECK_NEAR(b.getRotation().z, 0.0, 1e-3);
    }

    // D) A force applied off-centre creates torque in the right direction.
    {
        RigidBody b = makeBox(1, Vector3(0.5, 0.5, 0.5));
        b.calculateDerivedData();
        b.addForceAtPoint(Vector3(0, 1, 0), Vector3(1, 0, 0));    // τ = (1,0,0)×(0,1,0) = +z
        b.integrate(0.01);
        CHECK(b.getRotation().z > 0);
        CHECK_NEAR(b.getRotation().x, 0.0, 1e-9);
        CHECK_NEAR(b.getRotation().y, 0.0, 1e-9);
    }

    // E) World inverse inertia tensor after a 90° rotation about Z swaps x/y.
    {
        RigidBody b; b.setMass(1);
        Matrix3 I; I.setInertiaTensorCoeffs(2, 3, 4); b.setInertiaTensor(I);   // Iinv = diag(1/2,1/3,1/4)
        b.setOrientation(Quaternion(std::cos(real_pi / 4), 0, 0, std::sin(real_pi / 4)));
        b.calculateDerivedData();
        Matrix3 w = b.getInverseInertiaTensorWorld();
        CHECK_NEAR(w.data[0], 1.0 / 3, 1e-6);      // x now carries former-y inertia
        CHECK_NEAR(w.data[4], 1.0 / 2, 1e-6);
        CHECK_NEAR(w.data[8], 1.0 / 4, 1e-6);
        CHECK_NEAR(w.data[1], 0.0, 1e-6);
    }

    return test::report("rigidbody");
}
