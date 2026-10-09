#pragma once
// Stage 4: narrow-phase collision detection.
//
// collide(a, b) answers "do these two shapes overlap, and if so how do I describe the overlap?"
// The description is a Manifold: a contact normal plus one or two contact points.
//
// Conventions (everything later builds on these):
//   - the normal is a unit vector pointing from body A toward body B;
//   - moving B by normal * depth (or A by -normal * depth) just separates the pair;
//   - each contact point sits halfway through the overlap, between the two surfaces.
#include "body.hpp"

namespace phys {

struct ContactPoint {
    Vec2 point;   // world space
    Real depth;   // penetration at this point, measured along the normal (>= 0)
};

struct Manifold {
    Vec2 normal;       // from A to B
    Real depth = 0;    // minimum translation distance along `normal` that separates the shapes
    int count = 0;     // 1 or 2 when colliding
    ContactPoint points[2];
};

// Returns true and fills `out` when the bodies' shapes overlap. Works for any pairing of circle and
// convex polygon. `out` is untouched when it returns false.
bool collide(const Body& a, const Body& b, Manifold& out);

}  // namespace phys
