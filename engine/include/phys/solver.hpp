#pragma once
// Stage 5: contact response by sequential impulses (extended in stage 8 to solve joints alongside).
//
// Every contact point is a constraint on the relative velocity of the two bodies at that point:
//   normal:   v_n >= bias         (don't approach; bounce if we arrived fast; get out of any overlap)
//   tangent:  |impulse_t| <= mu * impulse_n   (Coulomb friction)
// We solve them one at a time, each by an impulse that fixes its own constraint given the current
// velocities. Fixing one disturbs its neighbours, so we sweep over all of them several times and the
// velocities converge. The key to making repeated sweeps correct is to clamp the ACCUMULATED impulse
// (sum over all sweeps), not each sweep's increment: a contact may push, but never pull.
//
// Joints (joints.hpp) are rows of the same kind, so they take part in the same sweeps: each sweep
// handles the joints first, then the contacts, since keeping a hinge together matters more than
// resolving a touch.
#include "body.hpp"
#include "collision.hpp"
#include "joints.hpp"
#include "solver_settings.hpp"

#include <vector>

namespace phys {

// One colliding pair found by the narrow phase: indices into the world's body array.
struct ContactPair {
    int a = 0, b = 0;
    Manifold manifold;
    // Both bodies are asleep (or static): the contact is remembered, with its impulses, so it can warm
    // start the pile when it wakes, but it takes no part in solving.
    bool dormant = false;
};

// Changes the velocities of the bodies so that no contact is moving into its partner (bouncing and
// applying friction as the materials dictate) and every joint's constraint holds. `dt` converts
// position error into a correcting speed. Reads each contact point's and joint row's stored impulses
// (when warm starting) and writes the final ones back.
void solve_constraints(std::vector<Body>& bodies, std::vector<ContactPair>& contacts, std::vector<Joint>& joints,
                       Real dt, const SolverSettings& settings);

// Copies solver impulses from last step's contacts to this step's, point by point: a point inherits
// the impulses of the previous point with the same body pair and the same feature id. Points with no
// match (a new touch) start from zero.
void transfer_impulses(const std::vector<ContactPair>& previous, std::vector<ContactPair>& current);

}  // namespace phys
