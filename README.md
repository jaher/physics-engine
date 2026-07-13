# A game physics engine (Millington's *Cyclone*, from scratch)

A complete C++ rigid-body physics engine implementing the architecture of Ian
Millington's **_Game Physics Engine Development_** — particles, mass-aggregate
systems, rigid bodies, coarse + fine collision detection, and the full
contact-resolution pipeline. Header-only, double precision, dependency-free, with
a physical test suite covering every layer.

![verification dashboard](docs/physics_demo.svg)

## What's implemented (mapped to the book)

| Part | Chapters | Feature | Files |
|------|----------|---------|-------|
| — | — | Vector3, Quaternion, Matrix3, Matrix4 (transforms, inverse, inertia) | `core.h`, `precision.h` |
| I | 3–4 | Particle + semi-implicit Euler integration | `particle.h` |
| II | 5–6 | Force generators: gravity, drag, spring, anchored spring, bungee, buoyancy, stiff (analytic) spring | `pfgen.h` |
| II | 7 | Particle contacts + resolver (velocity impulse & interpenetration, resting-contact correction) | `pcontacts.h` |
| II | 7 | Hard constraints: cables, rods, and anchored variants | `plinks.h` |
| II | 8 | Mass-aggregate particle world | `pworld.h` |
| III | 9–10 | Rigid body: 6-DOF integration, world inertia tensor, sleep system | `body.h` |
| III | 11 | Rigid-body force generators: gravity, spring, aero + control surface, buoyancy | `fgen.h` |
| IV | 12 | Broad phase: bounding-sphere hierarchy (BVH), potential contacts | `collide_coarse.h` |
| IV | 13 | Narrow phase: sphere/box/plane primitives, intersection tests, contact generation incl. **box-box SAT** | `collide_fine.h` |
| V | 14–15 | Contact + resolver: impulse velocity resolution, anisotropic **friction**, nonlinear position projection | `contacts.h` |
| V | 15 | Rigid-body world; ball-and-socket **joints** | `world.h`, `joints.h` |
| — | app. | Deterministic random source (vectors, quaternions) | `random.h` |
| + | — | **Cloth** — Verlet + structural/shear/bend constraints, sphere/ground collision with friction, wind | `cloth.h` |
| + | — | **Hair** — Verlet strands with segment + bend constraints, head collision, gusting wind + turbulence | `hair.h` |
| + | — | **Destruction** — grid/jitter fracture of a solid into welded box fragments that burst apart on impact | `fracture.h` |
| + | — | **Spring harmonics** — coupled spring-mass lattice (Hooke links via `ParticleSpring`, fixed ends) on the engine's particles with a symplectic step; excites the analytic normal modes / standing-wave harmonics `ω_n = 2√(k/m)·sin(nπ/2(P+1))` | `springs.h` |
| + | — | **Falling burning leaf** — tumbling-plate aerodynamics (anisotropic drag, reactive lift, underdamped broadside-seeking torque → flutter/tumble/spiral descent) coupled to a per-leaf combustion CA over an oak silhouette (ignite, glowing front, char, holes, curl) | `leaf.h` |
| + | — | **Non-Newtonian fluids** — weakly-compressible **SPH** (poly6/spiky/viscosity kernels, grid neighbour search, XSPH) with a *shear-rate-dependent* Herschel–Bulkley/power-law viscosity `μ(γ̇)` → shear-**thickening** (oobleck), shear-**thinning** (ketchup), **yield-stress**/Bingham (**lava**) and Newtonian fluids; two-way-coupled rigid balls, static cylinder obstacles, runtime emitters/drains | `sph.h` |
| + | — | **Screen-space fluid surface** — particles → eye-space depth → bilateral smoothing → normal reconstruction → shaded liquid (glossy Fresnel or molten blackbody emissive); makes SPH read as a surface, not spheres | `demos/common/fluidsurf.h` |
| + | — | **Framework parity** (vs Brax/Chrono/Gazebo/MuJoCo/ODE/PhysX/PyBullet/Webots/Unity — see [FEATURES.md](FEATURES.md)): capsules, raycasts, hinge/slider/fixed/distance/universal/gear joints + limits + motors + PD servos, heightfield terrain, CCD, contact events, triggers, kinematic bodies, RK4, raycast vehicle, character controller | `shapes.h`, `raycast.h`, `joints2.h`, `terrain.h`, `extras.h`, `vehicle.h`, `character.h` |

Include everything with `#include "phys/phys.h"`; the namespace is `phys`.  The
cloth and hair modules (extensions beyond the book, driving the 3D demos) are in
`phys/cloth.h` and `phys/hair.h`.

## Build & test

Header-only — just add `include/` to your include path. To build the tests and
the demo:

```sh
cmake -B build -S . && cmake --build build && ctest --test-dir build   # via CMake
bash tests/run_all.sh                                                  # or directly with g++
./build/demo docs/physics_demo.svg                                     # regenerate the dashboard
```

## Verification

Every layer is checked against analytic or known-physical results
(`tests/*.cpp`, **179 assertions, all passing**):

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
| `leaf` | the oak mask is a proper silhouette with veins; a dropped leaf falls at bounded speed, drifts sideways (flutters) and keeps an orthonormal frame; an unlit leaf keeps its fuel while an ignited one burns down to nothing |

The dashboard above is produced by `demos/demo.cpp` and shows, with no tuning:
decaying restitution bounces, quadratic drag shortening a projectile's range, a
tumbling box settling as its kinetic energy spikes-then-decays to zero, and a rod
pendulum tracing a constant-radius arc.

## 3D OpenGL demos

Three interactive real-time demos live in `demos/` and render with a small
modern-OpenGL pipeline (`demos/common/`): directional-light **shadow mapping**
(PCF), **Cook-Torrance PBR**, 4× **MSAA** and **ACES** tonemapping.

| Demo | What it shows |
|------|---------------|
| `rigid3d` | boxes and spheres dropping and piling up — the engine's rigid bodies, collision detection and contact resolver in 3D |
| `cloth3d` | a woven sheet draping over a sphere and pooling on the floor (`phys::Cloth`), two-sided fabric shading with sheen |
| `hair3d`  | ~3400 simulated guide strands → **24k rendered strands** via clumping (`phys::Hair`), **Kajiya-Kay** anisotropic shading, gusting wind |
| `clothesline3d` | garments pinned along a sagging string, hanging and swaying in a breeze (`phys::Cloth`) |
| `destruction3d` | a solid block shattered by a wrecking ball — **granite** chunks or **wood** splinters (`phys::fracture`) |
| `burn3d` | a hanging **antique world map** — a scan of the 1900 *Larousse* planisphere (public domain), texture-mapped onto the cloth — lit at one corner burns diagonally: a fire-propagation cellular automaton (`phys::BurningPaper`) drives a glowing ragged front that chars the print to black, with curling, holes and rising smoke + embers, while a light gusting breeze + travelling-wave `flutter` make the sheet undulate as it burns |
| `harmonics3d` | six coil springs strung between posts, each ringing in a pure normal mode (1…6 antinodes) of a coupled spring-mass chain (`phys::SpringChain`) — the **harmonic series** as standing waves, the higher modes visibly faster; each coil is a helix swept along the live masses ([`docs/spring_harmonics.mp4`](docs/spring_harmonics.mp4)) |
| `nonnewtonian3d` | the same heavy ball dropped into two SPH troughs — shear-thickening **oobleck** vs **water** (`phys::SPHFluid`): the oobleck's viscosity spikes under the impact's shear and absorbs the ball with barely a splash, while the water erupts in a crown. Drawn as a smooth screen-space surface ([`docs/nonnewtonian.mp4`](docs/nonnewtonian.mp4)) |
| `lava3d` | **hyper-real lava** — a hot shear-thinning Bingham fluid (`phys::Rheology::lava`) poured down a tilted channel and flowing around a stone **cylinder**, splitting and rejoining in its wake; a per-particle temperature drives a molten blackbody glow (yellow-white vent → dark crust) over the screen-space fluid surface ([`docs/lava.mp4`](docs/lava.mp4)) |
| `leaves3d` | **dozens of oak leaves spiralling down and burning** — tumbling-plate aerodynamics (`phys::FallingLeaf`) flutter each leaf down; it ignites partway and a combustion CA burns a glowing front across it, charring, curling and holing it before it is consumed, shedding smoke + embers. Leaves are **real oak-leaf photo textures** (a texture array; alpha defines each silhouette and drives the burn mask) ([`docs/burning_leaves.mp4`](docs/burning_leaves.mp4)) |
| `wall3d` | a metal ball punches **through** a concrete wall, releasing only the fragments inside a jagged angle-modulated radius — the wall survives with a ragged hole |
| `playground3d` | the framework-parity features in one scene: heightfield terrain, raycast **vehicle**, **character controller** hopping the dunes, motor-driven hinge **windmill**, **capsules** tumbling, and a **trigger zone** that lights up as the car passes ([`docs/playground.mp4`](docs/playground.mp4)) |

![rigid](docs/rigid3d.png)
![cloth](docs/cloth3d.png)
![hair](docs/hair3d.png)
![clothesline](docs/clothesline.png)
![granite](docs/destruction_granite.png)
![wood](docs/destruction_wood.png)

Rendered clips: [`cloth_drop.mp4`](docs/cloth_drop.mp4),
[`hair_wind.mp4`](docs/hair_wind.mp4), [`clothesline.mp4`](docs/clothesline.mp4),
[`destruction_granite.mp4`](docs/destruction_granite.mp4),
[`destruction_wood.mp4`](docs/destruction_wood.mp4),
[`wall_breakthrough.mp4`](docs/wall_breakthrough.mp4).

`destruction3d` takes `--material granite|wood`; the block is pre-split into
welded box fragments (`fractureBoxGrid`) that stay put until the ball's impact
calls `Destructible::shatter`, which bursts them apart with a radial + directional
impulse — from then on they are ordinary rigid bodies that collide and settle.

`burn3d` texture-maps a scan of the **1900 Larousse world planisphere**
("Planisphère des colonies européennes — Heure universelle", public domain, via
Wikimedia Commons) onto the burning cloth. The image is packed as zlib-compressed
RGB in `demos/assets/worldmap.pmz` and inflated at load, so run the demo from the
repo root (`./build/burn3d`) where that path resolves. The fire CA drives the char
(which blackens the print), holes and curl; the sheet texture is otherwise static.

`leaves3d` uses **real oak-leaf photo textures** (public domain / CC, via Wikimedia
Commons: a *Quercus rubra* autumn leaf plus *Quercus coccinea* scarlet-oak leaves
recoloured into autumn tints). They're cut out to alpha, packed as a zlib-compressed
RGBA texture array in `demos/assets/leaves.lefa`, and inflated at load — so run it
from the repo root. Each leaf's alpha both draws its silhouette and, sampled on a
grid, defines which cells can burn.

```sh
cmake --build build           # builds rigid3d, cloth3d, hair3d if GLFW/GLEW present
./build/rigid3d               # drag to orbit, scroll to zoom, Esc to quit
./build/cloth3d               # (hair3d likewise)
./build/cloth3d --shot out.png [frames]        # headless still
./build/cloth3d --video frames/f [nframes]     # headless image sequence for a video
```

Dependencies for the GL demos: **GLFW3, GLEW, GLM, zlib** (the engine and its
tests need none of these). Screenshots are written with a built-in zlib PNG
encoder; assemble a video from a `--video` sequence with any tool, e.g.
`ffmpeg -framerate 60 -i frames/f_%04d.png -pix_fmt yuv420p out.mp4`.

## Design notes / honest deviations

- **Header-only, `real = double`.** One `typedef` in `precision.h` switches to
  `float`; all maths goes through the `real_*` wrappers, as in the book.
- **Sleep system.** Bodies deactivate when their recency-weighted motion drops
  below `sleepEpsilon`. A body is constructed already above the threshold so it
  doesn't sleep on frame 1 (the book sets this when a body is woken).
- **Sphere-sphere contact point** uses the geometric overlap midpoint rather than
  the book's `positionOne + midline·½` (which can fall behind the first sphere) —
  more robust for the resolver's torque calculation.
- **Box-box** generates one contact (deepest feature) per pair, per the book, so
  large flat stacks rely on multiple frames to stabilise; point/edge and
  face/vertex cases are both handled via the 15-axis separating-axis test.
