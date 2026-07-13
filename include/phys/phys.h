// Umbrella header for the whole engine — a from-scratch C++ implementation of
// the "Cyclone" architecture from Ian Millington's *Game Physics Engine
// Development*, spanning particles, mass-aggregate systems, rigid bodies, coarse
// and fine collision detection, and the contact-resolution pipeline.
#pragma once

#include "precision.h"
#include "core.h"
#include "random.h"

// Particle physics (Parts I-II)
#include "particle.h"
#include "pfgen.h"
#include "pcontacts.h"
#include "plinks.h"
#include "pworld.h"

// Rigid-body physics (Part III)
#include "body.h"
#include "fgen.h"

// Collision detection (Part IV)
#include "collide_coarse.h"
#include "collide_fine.h"
#include "collide_convex.h"   // GJK/EPA convex hulls, cylinder & cone, static triangle mesh

// Contact physics (Part V)
#include "material.h"         // per-material friction/restitution + combine modes
#include "contacts.h"         // + persistent warm-started manifolds, rolling/spin friction
#include "joints.h"
#include "world.h"

// Destruction
#include "fracture.h"
#include "plastic.h"

// Oscillations
#include "springs.h"        // coupled spring-mass lattice → normal modes / harmonics
#include "leaf.h"           // falling burning leaf: tumbling-plate aero + combustion

// Fluids
#include "sph.h"            // non-Newtonian SPH (shear-thickening/thinning, yield stress)

// Feature parity with mainstream engines (Brax/Chrono/Gazebo/MuJoCo/ODE/PhysX/
// PyBullet/Webots/Unity): capsules, raycasts, joint zoo + motors, heightfields,
// CCD, contact events, triggers, kinematics, RK4, vehicles, characters.
#include "shapes.h"
#include "raycast.h"
#include "joints2.h"
#include "terrain.h"
#include "extras.h"
#include "vehicle.h"
#include "character.h"

// Reduced-coordinate dynamics, soft bodies, soft constraints & fields
#include "articulation.h"     // Featherstone O(n) articulated-body algorithm
#include "softbody.h"         // co-rotational tet-FEM soft bodies, PBF fluids, cloth self-collision
#include "constraint2.h"      // CFM/ERP soft link, breakable joint, conveyor, wind/radial fields, hydrodynamics, implicit spring

// Compound (multi-shape) rigid bodies
#include "compound.h"         // aggregate mass/COM/parallel-axis inertia over child primitives

// Scene queries, broad phases, serialization
#include "query.h"            // swept + overlap scene queries, speculative contacts
#include "broadphase2.h"      // sweep-and-prune + dynamic AABB tree
#include "serialize.h"        // rigid-body state snapshot / restore

// Robotics layer
#include "loader.h"           // URDF/MJCF model loader
#include "robotics.h"         // IMU/lidar/force sensors, RNE inverse dynamics, Jacobian IK, tendons

// Advanced / paradigm
#include "autodiff.h"         // forward-mode Dual → differentiable physics
#include "parallel.h"         // island-parallel solver (needs -pthread)
#include "decompose.h"        // approximate convex decomposition (V-HACD-lite)

// Electromagnetics
#include "fdtd.h"             // FDTD solution of Maxwell's curl equations on a Yee grid (TM/TE)
#include "mom.h"              // Method of Moments — thin-wire dipole EFIE (current, impedance, pattern)

// Middleware
#include "ros_bridge.h"       // ROS2-compatible messages + pub/sub graph (in-process; rclcpp-bindable)
// GPU execution lives in include/phys/gpu/gpu_sph.cuh — a CUDA module, compiled with nvcc (not this header-only path).
