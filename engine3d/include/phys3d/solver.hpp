#pragma once
// Stage 10, step 4: contact response in 3D by sequential impulses.
//
// The method is the 2D engine's (see phys/solver.hpp): every contact point is a constraint on the relative
// velocity of the two bodies there, each is solved by an impulse given the current velocities, all are
// swept several times, and the ACCUMULATED impulse is what gets clamped. Two things differ in 3D:
//   - friction acts in a plane, not along a line. Each point gets two tangent rows (t1, t2 perpendicular to
//     the normal and to each other), and it is the LENGTH of the combined tangent impulse that Coulomb's law
//     limits: |(jt1, jt2)| <= mu * jn, a circle. Clamping the two rows separately would make it a square,
//     with 41% more friction along the diagonals than along t1 or t2.
//   - how hard it is to accelerate a contact point involves the inertia TENSOR:
//         1/m_eff = 1/mA + 1/mB + dir . ((IA^-1 (ra x dir)) x ra) + dir . ((IB^-1 (rb x dir)) x rb).
#include "body.hpp"
#include "collision.hpp"

#include <phys/solver_settings.hpp>

#include <vector>

namespace phys3d {

using phys::SolverSettings;

// One colliding pair found by the narrow phase: indices into the world's body array.
struct ContactPair {
    int a = 0, b = 0;
    Manifold manifold;
    // Both bodies are asleep (or static): the contact is remembered, with its impulses, so it can warm
    // start the pile when it wakes, but it takes no part in solving.
    bool dormant = false;
};

// Changes the bodies' velocities so no contact is closing, with bounce and friction from their materials.
// Reads each point's stored impulses (when warm starting) and writes the final ones back.
void solve_contacts(std::vector<Body>& bodies, std::vector<ContactPair>& contacts, Real dt, const SolverSettings& settings);

// Copies impulses from last step's contacts to this step's: a point inherits those of the previous point
// with the same body pair and the same feature id. Unmatched points start from zero.
void transfer_impulses(const std::vector<ContactPair>& previous, std::vector<ContactPair>& current);

// Two unit vectors that, with the unit vector n, form a right-handed orthonormal frame (t1 x t2 = n). The
// choice depends only on n, so a contact whose normal holds steady keeps the same tangents from frame to
// frame, which is what lets its stored friction impulses be reused.
void tangent_basis(Vec3 n, Vec3& t1, Vec3& t2);

}  // namespace phys3d
