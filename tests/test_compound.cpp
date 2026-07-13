// Compound bodies: aggregate mass, centre of mass and the parallel-axis inertia
// tensor are correct, and the child primitives land at the right world offsets.
#include "phys/compound.h"
#include "phys/body.h"
#include "check.h"
#include <cmath>
using namespace phys;

int main() {
    // A) symmetric barbell: two spheres (m=1, r=0.1) at x=±0.5
    {
        RigidBody b; CompoundShape cs;
        cs.addSphere(0.1, Vector3(-0.5, 0, 0), 1.0);
        cs.addSphere(0.1, Vector3(0.5, 0, 0), 1.0);
        cs.finalize(&b);
        CHECK_NEAR(cs.totalMass, 2.0, 1e-9);
        CHECK_NEAR(cs.com.x, 0.0, 1e-9);
        // I along the bar axis (x) = 2·(2/5·m·r²);  ⟂ axis (z) also gets 2·m·d²
        CHECK_NEAR(cs.inertia.data[0], 2 * (0.4 * 1 * 0.01), 1e-9);
        CHECK_NEAR(cs.inertia.data[8], 2 * (0.4 * 1 * 0.01 + 1 * 0.25), 1e-9);
        CHECK(cs.spheres.size() == 2);
    }

    // B) asymmetric masses shift the COM
    {
        RigidBody b; CompoundShape cs;
        cs.addSphere(0.1, Vector3(-0.5, 0, 0), 1.0);
        cs.addSphere(0.1, Vector3(0.5, 0, 0), 3.0);
        cs.finalize(&b);
        CHECK_NEAR(cs.totalMass, 4.0, 1e-9);
        CHECK_NEAR(cs.com.x, (1 * -0.5 + 3 * 0.5) / 4.0, 1e-9);       // = 0.25
    }

    // C) child primitives are placed relative to the COM in the body frame
    {
        RigidBody b; CompoundShape cs;
        cs.addBox(Vector3(0.2, 0.2, 0.2), Vector3(0, 0, 0), 1.0);
        cs.addBox(Vector3(0.2, 0.2, 0.2), Vector3(0, 0.4, 0), 1.0);   // stacked → COM at y=0.2
        cs.finalize(&b);
        CHECK_NEAR(cs.com.y, 0.2, 1e-9);
        b.setPosition(Vector3(0, 0, 0)); b.setOrientation(Quaternion()); b.calculateDerivedData();
        cs.calculateInternals();
        // lower box sits at y = 0 - 0.2 (below COM), upper at +0.2, in world (body at origin)
        CHECK_NEAR(cs.boxes[0]->getAxis(3).y, -0.2, 1e-6);
        CHECK_NEAR(cs.boxes[1]->getAxis(3).y, 0.2, 1e-6);
        // combined inertia about COM exceeds a single box's (mass is spread out)
        Matrix3 single; single.setBlockInertiaTensor(Vector3(0.2, 0.2, 0.2), 1.0);
        CHECK(cs.inertia.data[0] > single.data[0]);
    }

    // D) a compound body integrates as one rigid body (falls under gravity)
    {
        RigidBody b; CompoundShape cs;
        cs.addBox(Vector3(0.2, 0.1, 0.1), Vector3(-0.3, 0, 0), 1.0);
        cs.addBox(Vector3(0.2, 0.1, 0.1), Vector3(0.3, 0, 0), 1.0);
        cs.finalize(&b);
        b.setPosition(Vector3(0, 5, 0)); b.setAcceleration(0, -9.81, 0); b.setDamping(1, 1); b.calculateDerivedData();
        for (int i = 0; i < 100; i++) { b.clearAccumulators(); b.integrate(0.01); }   // 1 s free fall
        CHECK_NEAR(b.getPosition().y, 5.0 - 0.5 * 9.81 * 1.0, 0.2);
    }

    return test::report("compound");
}
