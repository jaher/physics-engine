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

// Contact physics (Part V)
#include "contacts.h"
#include "joints.h"
#include "world.h"

// Destruction
#include "fracture.h"
#include "plastic.h"

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
