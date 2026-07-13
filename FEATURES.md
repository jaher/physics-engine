# Framework parity matrix

Feature survey of nine mainstream physics frameworks — **Brax, Chrono, Gazebo,
MuJoCo, ODE, PhysX, PyBullet, Webots, Unity** — and this engine's status on each
feature family. Everything marked ✅ new is compile-tested and covered by
`tests/test_features.cpp` (52 assertions).

| Feature family | Found in | Status here | Where |
|---|---|---|---|
| Rigid bodies (6-DOF, inertia tensors) | all nine | ✅ (core) | `body.h` |
| Sphere / box / plane primitives | all nine | ✅ (core) | `collide_fine.h` |
| **Capsule primitive** + all pair tests | all nine | ✅ **new** | `shapes.h` |
| **Raycasting** (sphere/box/plane/capsule, world query) | all nine | ✅ **new** | `raycast.h` |
| **Hinge/revolute joint** + limits + **velocity motor** | ODE, MuJoCo, PhysX, Bullet, Chrono, Unity, Gazebo, Webots | ✅ **new** | `joints2.h` |
| **Slider/prismatic joint** + limits + motor | same | ✅ **new** | `joints2.h` |
| **Fixed / distance / universal joints** | same | ✅ **new** | `joints2.h` |
| **Gear coupling** | ODE, Chrono | ✅ **new** | `joints2.h` |
| **PD position servo** (actuator) | MuJoCo, Webots, Gazebo | ✅ **new** (`HingeServo`) | `joints2.h` |
| Ball-and-socket joint | all | ✅ (core + new PBD version) | `joints.h`, `joints2.h` |
| **Heightfield terrain** (+ raycast) | ODE, PhysX, Bullet, Gazebo, Unity | ✅ **new** | `terrain.h` |
| **Continuous collision detection** (swept sphere) | PhysX, Bullet, Unity, Chrono | ✅ **new** | `extras.h` |
| **Contact events** (begin/end callbacks) | Unity, PhysX, Gazebo, Webots | ✅ **new** | `extras.h` |
| **Trigger volumes** | Unity, PhysX | ✅ **new** | `extras.h` |
| **Kinematic (animated) bodies** | all nine | ✅ **new** helpers | `extras.h` |
| **RK4 integrator** | Brax, Chrono | ✅ **new** | `extras.h` |
| **Raycast vehicle** (suspension, engine, steering, tyre grip) | Bullet, PhysX, Unity, Chrono, Webots | ✅ **new** | `vehicle.h` |
| **Character controller** (capsule collide-and-slide, jump, grounded) | PhysX, Unity, Bullet | ✅ **new** | `character.h` |
| Broad phase (BVH) | all nine | ✅ (core) | `collide_coarse.h` |
| Sleeping / deactivation | all nine | ✅ (core) | `body.h` |
| Springs, drag, buoyancy, aero | most | ✅ (core) | `fgen.h`, `pfgen.h` |
| Cloth | PhysX, Bullet, MuJoCo(flex) | ✅ | `cloth.h` |
| Fracture / destruction | PhysX (blast), Chrono | ✅ | `fracture.h` |
| Determinism | Brax, MuJoCo, ODE | ✅ (fixed-step, no hidden state) | — |

**Not implemented** (noted honestly): Featherstone articulated-body dynamics
(MuJoCo/PhysX reduced coordinates — our joints are maximal-coordinate iterative),
convex-hull & triangle-mesh colliders, soft-body volumes/FEM (Chrono), fluids,
GPU/autodiff simulation (Brax's differentiability), and the robotics middleware
layers (ROS integration, sensors, SDF/URDF loading) of Gazebo/Webots — those are
toolchain features rather than physics.

Demo: `demos/playground3d.cpp` exercises the new features in one scene
(heightfield + vehicle + character + motorized hinge windmill + capsules +
trigger zone) — see `docs/playground.mp4`.
