#pragma once
// Stage 10, step 3: narrow-phase collision detection in 3D.
//
// Same contract as the 2D engine's collide():
//   - the normal is a unit vector pointing from body A toward body B;
//   - moving B by normal * depth (or A the other way) separates the pair;
//   - each contact point sits halfway through the overlap, between the two surfaces.
// What changes in 3D is the size of a manifold. Two faces meeting share a polygon, not a segment, so a
// manifold holds up to FOUR points: three are needed to stop a box rocking on a face, and four keep the
// support polygon as large as possible. And two boxes can touch edge against edge, a case with no 2D
// counterpart, which gives a single point.
#include "body.hpp"

#include <cstdint>

namespace phys3d {

struct ContactPoint {
    Vec3 point;      // world space
    Real depth = 0;  // penetration at this point along the normal (>= 0)
    // Names the pair of features (faces, edges, vertices) that produced the point, so it can be matched
    // with the same point in the previous frame. Only comparable within one body pair.
    std::uint32_t id = 0;
    // Solver state carried across frames (filled in by the solver, never by collide()).
    Real normal_impulse = 0;
    Real tangent_impulse[2] = {0, 0};
};

struct Manifold {
    static constexpr int kMaxPoints = 4;

    Vec3 normal;     // from A to B
    Real depth = 0;  // how far B must move along `normal` to separate the pair
    int count = 0;   // 1 to 4 when colliding
    ContactPoint points[kMaxPoints];
};

// True and fills `out` when the bodies' shapes overlap; `out` is untouched otherwise. Handles every
// pairing of sphere and box.
bool collide(const Body& a, const Body& b, Manifold& out);

// Is `world_point` inside the body's shape? (For picking, and for testing the narrow phase.)
bool contains(const Body& body, Vec3 world_point);

}  // namespace phys3d
