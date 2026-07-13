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

## Second wave — closing the gap matrix

A follow-up audit (`docs` gap analysis) drove a second round that implemented
essentially every remaining capability family. All are header-only and covered by
their own test suites (`bash tests/run_all.sh` → **434 assertions across 18 suites**).

| Was missing in | Now here | Where | Test |
|---|---|---|---|
| Convex-hull & triangle-mesh colliders, cylinder/cone (all nine) | ✅ GJK + EPA convex-convex, cyl/cone vs plane, sphere/box vs trimesh | `collide_convex.h` | `convex` (36) |
| Featherstone reduced-coordinate articulation (MuJoCo/PhysX/Bullet/Brax) | ✅ O(n) Articulated-Body Algorithm, revolute/prismatic serial-tree chains | `articulation.h` | `articulation` (6) |
| Persistent multi-point manifolds + warm starting (Bullet/PhysX) | ✅ up-to-4-point cache, warm-started impulses, material combine modes, rolling/spin friction | `contacts.h`, `material.h` | `contacts2` (17) |
| Volumetric FEM soft bodies (Chrono/PhysX5), cloth self-collision, PBF fluids | ✅ co-rotational tet-FEM, spatial-hash self-collision, position-based fluids | `softbody.h` | `softbody` (22) |
| Soft constraints (CFM/ERP), breakable joints, conveyors, force fields, hydrodynamics, implicit integrator | ✅ all | `constraint2.h` | `constraint2` (45) |
| Shape sweeps & overlap queries, SAP / dynamic-AABB-tree broad phases, serialization, speculative CCD | ✅ all | `query.h`, `broadphase2.h`, `serialize.h` | `query` (49) |
| URDF/MJCF loading, IMU/lidar/force sensors, inverse dynamics, IK, tendons (Gazebo/Webots/MuJoCo) | ✅ all | `loader.h`, `robotics.h` | `robotics` (53) |
| Differentiable physics (Brax/MJX), multithreaded solver, convex decomposition | ✅ forward-mode autodiff, island-parallel `std::thread` solver, V-HACD-lite | `autodiff.h`, `parallel.h`, `decompose.h` | `advanced` (27) |

**Still not implemented** (need external toolchains, not header-only physics):
**GPU execution** of the simulation (CUDA / a compute backend — the autodiff path
is differentiable but runs on the CPU) and the **ROS/ROS2 middleware** bridge
(the sensors, loaders and inverse dynamics exist; the pub/sub transport does not).

Demo: `demos/playground3d.cpp` exercises the new features in one scene
(heightfield + vehicle + character + motorized hinge windmill + capsules +
trigger zone) — see `docs/playground.mp4`.
