# A from-scratch C++ physics & simulation engine

A comprehensive, **header-only** C++ physics engine — written from scratch, double
precision, with a **dependency-free core** — that reaches well beyond rigid bodies.
It spans the full rigid-body dynamics and contact-resolution pipeline; convex &
triangle-mesh collision via **GJK + EPA**; joints, motors and constraints; **cloth
and hair**; **fracture, destruction and explosions**; **SPH fluids** — from a CPU
non-Newtonian rheology solver to a **CUDA GPU** solver that scales to **millions of
particles**; **FEM soft bodies**; **Featherstone articulated bodies**; a **robotics**
layer (URDF/MJCF, inverse dynamics, IK, sensors); differentiable and island-parallel
solvers; and even a computational-**electromagnetics** layer (FDTD + Method of
Moments).

Every layer is validated against analytic or known-physical results — **589
assertions across 23 test suites** (including a CUDA GPU suite) — and showcased by
**40+ real-time OpenGL demos** with shadow-mapped Cook-Torrance PBR rendering.

![verification dashboard](docs/physics_demo.svg)

## Gallery

Each thumbnail links to a rendered clip.

### Fluids — CPU rheology & CUDA GPU (to millions of particles)

<table>
<tr>
<td align="center" width="33%"><a href="docs/waterfall.mp4"><img src="docs/waterfall.png" width="270"></a><br><b>waterfall3d</b><br><sub>2.4M-particle GPU SPH waterfall</sub></td>
<td align="center" width="33%"><a href="docs/coriolis.mp4"><img src="docs/coriolis.png" width="270"></a><br><b>coriolis3d</b><br><sub>Coriolis effect: rotating vs still tank</sub></td>
<td align="center" width="33%"><a href="docs/gpu_fluid.mp4"><img src="docs/gpu_fluid.png" width="270"></a><br><b>gpu_fluid3d</b><br><sub>GPU dam-break, ~900k particles</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="docs/lava.mp4"><img src="docs/lava.png" width="270"></a><br><b>lava3d</b><br><sub>hot shear-thinning Bingham fluid</sub></td>
<td align="center" width="33%"><a href="docs/nonnewtonian.mp4"><img src="docs/nonnewtonian.png" width="270"></a><br><b>nonnewtonian3d</b><br><sub>shear-thickening oobleck vs water</sub></td>
<td align="center" width="33%"><a href="docs/snow_footprints.mp4"><img src="docs/snow_footprints.png" width="270"></a><br><b>snow3d</b><br><sub>footprints pressed into a drift</sub></td>
</tr>
</table>

### Fracture, destruction & explosions

<table>
<tr>
<td align="center" width="33%"><a href="docs/explosion.mp4"><img src="docs/explosion.png" width="270"></a><br><b>explosion3d</b><br><sub>Voronoi-fractured block + fireball</sub></td>
<td align="center" width="33%"><a href="docs/bulletbarrel.mp4"><img src="docs/bulletbarrel.png" width="270"></a><br><b>bulletbarrel3d</b><br><sub>bullet → barrel chain reaction</sub></td>
<td align="center" width="33%"><a href="docs/roomblast.mp4"><img src="docs/roomblast.png" width="270"></a><br><b>roomblast3d</b><br><sub>a room of objects blown apart</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="docs/bulletglass.mp4"><img src="docs/bulletglass.png" width="270"></a><br><b>bulletglass3d</b><br><sub>bullet through a shattering pane</sub></td>
<td align="center" width="33%"><a href="docs/bulletproofglass.mp4"><img src="docs/bulletproofglass.png" width="270"></a><br><b>bulletproofglass3d</b><br><sub>laminated glass stops the round</sub></td>
<td align="center" width="33%"><a href="docs/mirror_shatter.mp4"><img src="docs/mirror_shatter.png" width="270"></a><br><b>mirror3d</b><br><sub>a shattering mirror</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="docs/destruction_granite.mp4"><img src="docs/destruction_granite.png" width="270"></a><br><b>destruction3d</b><br><sub>granite chunks off a wrecking ball</sub></td>
<td align="center" width="33%"><a href="docs/destruction_wood.mp4"><img src="docs/destruction_wood.png" width="270"></a><br><b>destruction3d</b><br><sub>wood splinters</sub></td>
<td align="center" width="33%"><a href="docs/wall_breakthrough.mp4"><img src="docs/wall_breakthrough.png" width="270"></a><br><b>wall3d</b><br><sub>a ball punches a ragged hole</sub></td>
</tr>
</table>

### Rigid bodies, cloth, hair & soft bodies

<table>
<tr>
<td align="center" width="33%"><a href="docs/stacking.mp4"><img src="docs/stacking.png" width="270"></a><br><b>stack3d</b><br><sub>warm-started brick pyramid</sub></td>
<td align="center" width="33%"><a href="docs/convex.mp4"><img src="docs/convex.png" width="270"></a><br><b>convex3d</b><br><sub>convex-hull gems (GJK/EPA)</sub></td>
<td align="center" width="33%"><img src="docs/rigid3d.png" width="270"><br><b>rigid3d</b><br><sub>boxes & spheres piling up</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="docs/cloth_drop.mp4"><img src="docs/cloth3d.png" width="270"></a><br><b>cloth3d</b><br><sub>a sheet draping over a sphere</sub></td>
<td align="center" width="33%"><a href="docs/hair_wind.mp4"><img src="docs/hair3d.png" width="270"></a><br><b>hair3d</b><br><sub>24k rendered strands, Kajiya-Kay</sub></td>
<td align="center" width="33%"><a href="docs/clothesline.mp4"><img src="docs/clothesline.png" width="270"></a><br><b>clothesline3d</b><br><sub>garments swaying in a breeze</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="docs/softbody.mp4"><img src="docs/softbody.png" width="270"></a><br><b>softbody3d</b><br><sub>co-rotational FEM jelly cubes</sub></td>
<td align="center" width="33%"><a href="docs/rope.mp4"><img src="docs/rope.png" width="270"></a><br><b>rope3d</b><br><sub>rope dynamics</sub></td>
<td align="center" width="33%"><a href="docs/tangle.mp4"><img src="docs/tangle.png" width="270"></a><br><b>tangle3d</b><br><sub>tangled strands</sub></td>
</tr>
</table>

### Mechanisms, articulated bodies & robotics

<table>
<tr>
<td align="center" width="33%"><a href="docs/articulation.mp4"><img src="docs/articulation.png" width="270"></a><br><b>articulation3d</b><br><sub>Featherstone multi-pendulum</sub></td>
<td align="center" width="33%"><a href="docs/robotics.mp4"><img src="docs/robotics.png" width="270"></a><br><b>robotics3d</b><br><sub>arm IK tracking + lidar scan</sub></td>
<td align="center" width="33%"><a href="docs/playground.mp4"><img src="docs/playground.png" width="270"></a><br><b>playground3d</b><br><sub>vehicle + character + terrain</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="docs/galton.mp4"><img src="docs/galton.png" width="270"></a><br><b>galton3d</b><br><sub>Galton board → bell curve</sub></td>
<td align="center" width="33%"><a href="docs/constraints.mp4"><img src="docs/constraints.png" width="270"></a><br><b>constraints3d</b><br><sub>conveyor belt + crosswind field</sub></td>
<td align="center" width="33%"><a href="docs/spring_harmonics.mp4"><img src="docs/spring_harmonics.png" width="270"></a><br><b>harmonics3d</b><br><sub>standing-wave spring modes</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="docs/chains.mp4"><img src="docs/chains.png" width="270"></a><br><b>chains3d</b><br><sub>swinging chains</sub></td>
<td align="center" width="33%"><a href="docs/forge_bend.mp4"><img src="docs/forge_bend.png" width="270"></a><br><b>forge3d</b><br><sub>bending hot metal</sub></td>
<td align="center" width="33%"></td>
</tr>
</table>

### Combustion, aerodynamics & electromagnetics

<table>
<tr>
<td align="center" width="33%"><a href="docs/burning_paper.mp4"><img src="docs/burning_paper.png" width="270"></a><br><b>burn3d</b><br><sub>an antique world map burning</sub></td>
<td align="center" width="33%"><a href="docs/burning_leaves.mp4"><img src="docs/burning_leaves.png" width="270"></a><br><b>leaves3d</b><br><sub>oak leaves fluttering down & burning</sub></td>
<td align="center" width="33%"><a href="docs/antenna3d.mp4"><img src="docs/antenna3d.png" width="270"></a><br><b>antenna3d</b><br><sub>3-D dipole radiation pattern (MoM)</sub></td>
</tr>
<tr>
<td align="center" width="33%"><a href="docs/em_antenna.mp4"><img src="docs/em_antenna.png" width="270"></a><br><b>em2d</b><br><sub>Maxwell FDTD driven by a MoM dipole</sub></td>
<td align="center" width="33%"></td>
<td align="center" width="33%"></td>
</tr>
</table>

## Capabilities

Everything is header-only; include it all with `#include "phys/phys.h"` (namespace
`phys`). See [FEATURES.md](FEATURES.md) for a feature-by-feature comparison against
nine mainstream frameworks (Brax, Chrono, Gazebo, MuJoCo, ODE, PhysX, PyBullet,
Webots, Unity).

| Capability | Files |
|---|---|
| **Core maths** — `Vector3`, `Quaternion`, `Matrix3`, `Matrix4` (transforms, inverse, inertia tensors), deterministic random source | `core.h`, `precision.h`, `random.h` |
| **Particles & mass-aggregate** — semi-implicit Euler particles; force generators (gravity, drag, spring, anchored/bungee spring, buoyancy, stiff analytic spring); particle contacts + resolver; cables & rods; mass-aggregate world | `particle.h`, `pfgen.h`, `pcontacts.h`, `plinks.h`, `pworld.h` |
| **Rigid-body dynamics** — 6-DOF integration, world inertia tensor, sleep/deactivation system; rigid-body force generators (gravity, spring, aero + control surface, buoyancy) | `body.h`, `fgen.h` |
| **Collision detection** — bounding-sphere hierarchy (BVH) broad phase; sphere/box/plane narrow phase incl. **box-box SAT**; impulse contact resolver with anisotropic **friction** and nonlinear position projection | `collide_coarse.h`, `collide_fine.h`, `contacts.h` |
| **Convex & mesh geometry** — arbitrary convex hulls (**GJK** boolean + closest distance + **EPA** penetration), analytic cylinder & cone, static triangle-mesh collider with grid midphase | `collide_convex.h` |
| **Joints, motors & constraints** — ball-and-socket, hinge/slider/fixed/distance/universal, gear coupling, limits + velocity motors + PD servos; CFM/ERP soft links, breakable joints, conveyors, wind/radial force fields, rigid hydrodynamics, implicit spring | `joints.h`, `joints2.h`, `constraint2.h` |
| **Contact solver upgrade** — persistent up-to-4-point manifolds, warm-started impulses, `Material` combine modes, rolling & spinning friction | `contacts.h`, `material.h` |
| **Articulated bodies** — Featherstone O(n) articulated-body algorithm in reduced (joint) coordinates for serial/tree chains of revolute/prismatic joints | `articulation.h` |
| **Cloth & hair** — Verlet cloth (structural/shear/bend constraints, collision + friction, wind); Verlet hair strands with bend constraints, head collision, gusting turbulence | `cloth.h`, `hair.h` |
| **Soft bodies & PBF** — co-rotational tetrahedral **FEM** deformables, cloth self-collision, position-based fluids | `softbody.h` |
| **Fracture, destruction & explosions** — grid/jitter box fracture, **Voronoi** irregular convex chunks, radial+concentric **impact ("bullet-hole") glass** fracture, **hollow-cylinder shell** fracture (regular grid or Voronoi-irregular metal shards) — all with exact volume/COM/inertia; `detonate()` a charge for a radial blast with an upward plume | `fracture.h`, `voronoi.h`, `glassfrac.h`, `shellfrac.h` |
| **Non-Newtonian SPH fluids** — weakly-compressible SPH (poly6/spiky/viscosity kernels, grid neighbours, XSPH) with a shear-rate-dependent Herschel–Bulkley/power-law viscosity → shear-**thickening** (oobleck), shear-**thinning** (ketchup), **yield-stress** (lava) and Newtonian fluids; two-way-coupled rigid balls, obstacles, runtime emitters/drains | `sph.h` |
| **GPU fluids (CUDA)** — massively-parallel SPH (spatial-hash grid, thrust sort, on-GPU density/force/integrate) with **surface-tension cohesion**, **box (SDF) obstacles**, a **recirculating emitter**, a rotating-frame **Coriolis** force and advected **passive dye**; ~360k particles at ~0.2 ms/step on an RTX 5090, driving a 2.4M-particle waterfall | `gpu/gpu_sph.cuh` |
| **Screen-space fluid surface** — particle depth → bilateral smoothing → normal reconstruction → shaded liquid (glossy Fresnel or molten blackbody emissive) | `demos/common/fluidsurf.h` |
| **Aerodynamics & combustion** — tumbling-plate falling-leaf aerodynamics + per-leaf combustion CA; paper/cloth fire propagation with char, holes and smoke | `leaf.h`, and demo `burn3d` |
| **Spring harmonics** — coupled spring-mass lattice exciting analytic normal modes `ω_n = 2√(k/m)·sin(nπ/2(P+1))` | `springs.h` |
| **Robotics layer** — URDF/MJCF loader, IMU/lidar/contact-force sensors, recursive-Newton-Euler inverse dynamics, Jacobian IK, tendon/muscle actuators | `loader.h`, `robotics.h` |
| **Framework-parity extras** — capsules, raycasts, heightfield terrain, continuous collision detection, contact events, trigger volumes, kinematic bodies, RK4, raycast vehicle, character controller | `shapes.h`, `raycast.h`, `terrain.h`, `extras.h`, `vehicle.h`, `character.h` |
| **Scene queries & broad phases** — swept + overlap queries, speculative contacts, sweep-and-prune, dynamic AABB tree, state snapshot/restore | `query.h`, `broadphase2.h`, `serialize.h` |
| **Advanced paradigms** — forward-mode autodiff → differentiable physics, island-parallel `std::thread` solver, approximate convex decomposition (V-HACD-lite) | `autodiff.h`, `parallel.h`, `decompose.h` |
| **Compound (multi-shape) bodies** — one rigid body over several offset box/sphere primitives, with aggregate mass, COM and parallel-axis inertia | `compound.h` |
| **ROS/ROS2 bridge** — ROS2-compatible messages + node/publisher/subscription graph over an in-process transport, bindable to real `rclcpp` | `ros_bridge.h` |
| **Computational electromagnetics** — Maxwell's equations by Yee-grid **FDTD** (TM & TE, PEC/absorbing borders); dipole antenna analysis by the **Method of Moments** (thin-wire EFIE) | `fdtd.h`, `mom.h` |

## Build & run

Header-only — just put `include/` on your include path. To build the tests and demos:

```sh
cmake -B build -S . && cmake --build build && ctest --test-dir build   # via CMake
bash tests/run_all.sh                                                  # or directly with g++
```

The **3D OpenGL demos** need **GLFW3, GLEW, GLM, zlib** (the engine and its tests
need none of these); the **GPU** module needs **CUDA / nvcc**. Both are optional and
skipped by CMake if absent. Each demo is interactive (drag to orbit, scroll to zoom)
and also runs headless to render stills or image sequences:

```sh
./build/rigid3d                                # interactive
./build/waterfall3d --shot out.png [frames]    # headless still
./build/waterfall3d --video frames/f [nframes] # headless image sequence
```

Screenshots are written with a built-in zlib PNG encoder; assemble a video from a
`--video` sequence with any tool, e.g. `ffmpeg -framerate 30 -i frames/f_%04d.png
-pix_fmt yuv420p out.mp4`. `burn3d` and `leaves3d` load packed public-domain image
assets (an 1900 *Larousse* world planisphere; real oak-leaf photos) from
`demos/assets/`, so run them from the repo root.

## Verification

Every layer is checked against analytic or known-physical results (`tests/*.cpp`,
**589 assertions, all passing** — 23 suites, incl. a CUDA GPU suite):

| Suite | What it proves |
|-------|----------------|
| `core` | dot/cross, matrix inverse `M·M⁻¹=I`, quaternion→matrix rotation, transform round-trips, quaternion-integration convergence |
| `particle` | projectile matches `v=v₀+at`; spring settles at `rest+mg/k`; rod holds length; ball rebounds to `e²h`; cable never over-extends |
| `rigidbody` | free fall; torque gives `I⁻¹τ`; torque-free spin is conserved; off-centre force makes the right torque; world inertia tensor rotates correctly |
| `collision` | sphere/sphere, sphere/plane, box/plane (4 corners), box/sphere, **box/box** penetration & normals; BVH finds the overlapping pair only |
| `resolution` | a box settles at the right height; bounces; friction decelerates a slide; a sphere rests; an elastic bounce conserves height |
| `stacking` | interpenetrating boxes separate; a box comes to rest on a fixed box (box-box in the pipeline) |
| `harmonics` | chain mode frequencies match `ω_n = 2√(k/m)·sin(nπ/2(P+1))`; a single mass gives `√(2k/m)`; modes are ordered; an undamped chain conserves energy; a tensioned chain gives a near-linear overtone series |
| `sph` | the `μ(γ̇)` rheology thickens/thins/yields correctly; an SPH pool settles, conserves particle count and stays inside its box with no NaNs; a shear-thickening pool holds a dropped ball higher than water; runtime add/remove stays in sync and a cylinder ejects particles from its footprint |
| `convex` | GJK boolean + distance and EPA penetration on convex hulls (aligned/offset/rotated); sphere & box rest on a triangle-mesh ramp; cylinder/cone vs plane |
| `articulation` | a Featherstone pendulum's period = `2π√(I/mgl)`; released-link acceleration = `mgl/I`; a free 2-link chain conserves energy; a holding torque keeps a link static |
| `contacts2` | material combine modes; a warm-started stack settles with less penetration in fewer iterations; rolling friction brings a rolling sphere to rest |
| `softbody` | polar decomposition recovers rotation/identity; a dropped soft box squashes and rebounds keeping its volume; self-collision separates overlapping particles; PBF holds rest density |
| `constraint2` | CFM softness lowers restoring force; a breakable joint fires exactly at its load; a conveyor drives a box to belt speed; hydrodynamic drag gives a bounded terminal velocity; the implicit spring stays stable where explicit diverges |
| `query` | a swept sphere catches a target it would tunnel past; sweep-and-prune and the dynamic AABB tree match brute-force pairs; snapshot→restore is bitwise-exact |
| `robotics` | URDF/MJCF parse to the right links/joints; IMU reads `g` when static; RNE gravity torque = `mgl`; Jacobian IK error decreases monotonically; lidar returns the expected ranges |
| `advanced` | autodiff gradients match finite differences; the island-parallel solver is bit-identical to serial; an L-shape decomposes into convex pieces covering its area |
| `gpu` (CUDA) | the GPU grid density equals the CPU brute-force reference; a dropped block stays finite/in-bounds near rest density; scales to >200k particles at <60 ms/step; a solid **box obstacle** is never penetrated; the **recirculating emitter** conserves particle count and leaves none below the drain plane; **surface-tension cohesion** makes a free blob more compact; in a rotating frame the **Coriolis** force conserves a particle's speed and rotates its velocity at `2·|Ω|`, while an inertial-frame particle flies straight |
| `ros` | the node graph delivers every published message to all subscribers; Imu/LaserScan/JointState round-trip; sim→ROS adapters build well-formed messages |
| `compound` | a multi-shape body aggregates total mass, centre of mass and the parallel-axis inertia tensor (barbell/stack); child primitives land at the right COM-relative world offsets; it integrates as one rigid body |
| `em` | FDTD reproduces a PEC-cavity eigenfrequency to within 0.01% of `(c/2)√((1/Lx)²+(1/Ly)²)` and a wavefront travels at 0.9975 c; the leapfrog stays energy-bounded for 3000 steps; the MoM dipole gives a near-resonant Rin ≈ 73 Ω, a symmetric feed-peaked current, capacitive→inductive reactance, and a radiation pattern with a deep axial null |
| `explosion` | grid fracture conserves the block's volume and fragment count exactly; **Voronoi fracture** partitions the block into convex cells with valid mass properties; the radial+concentric **glass fracture** tiles the pane (minus the punched hole) into thin convex shards, and the **hollow-cylinder shell fracture** dices the barrel wall (grid or **irregular Voronoi**) into convex metal pieces that tile the shell; `detonate()` throws every fragment radially outward with a net upward plume momentum, its blast energy scales with the square of the charge strength, and fragments outside the blast radius stay asleep |
| `leaf` | the oak mask is a proper silhouette with veins; a dropped leaf falls at bounded speed, drifts sideways (flutters) and keeps an orthonormal frame; an unlit leaf keeps its fuel while an ignited one burns down to nothing |

The verification dashboard above is produced by `demos/demo.cpp` (`./build/demo
docs/physics_demo.svg`) and shows, with no tuning: decaying restitution bounces,
quadratic drag shortening a projectile's range, a tumbling box settling as its
kinetic energy spikes-then-decays to zero, and a rod pendulum tracing a
constant-radius arc.

## Rendering

The 3D demos share a small modern-OpenGL pipeline (`demos/common/`): directional-light
**shadow mapping** (PCF), **Cook-Torrance PBR**, 4× **MSAA** and **ACES** tonemapping,
plus the screen-space fluid surface for SPH liquids.

## Design notes

- **Header-only, `real = double`.** One `typedef` in `precision.h` switches the whole
  engine to `float`; all maths goes through the `real_*` wrappers.
- **Sleep system.** Bodies deactivate when their recency-weighted motion drops below
  `sleepEpsilon`, and are constructed already above the threshold so they don't sleep
  on frame 1.
- **Sphere-sphere contact point** uses the geometric overlap midpoint (robust for the
  resolver's torque calculation).
- **Box-box** generates one contact (deepest feature) per pair via a 15-axis
  separating-axis test, so large flat stacks stabilise over a few frames; the
  warm-started manifold solver (`contacts2`) is used where persistent multi-point
  contact matters.
</content>
