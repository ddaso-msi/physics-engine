#pragma once
// Stage 5: contact response by sequential impulses.
//
// Every contact point is a constraint on the relative velocity of the two bodies at that point:
//   normal:   v_n >= bias         (don't approach; bounce if we arrived fast; get out of any overlap)
//   tangent:  |impulse_t| <= mu * impulse_n   (Coulomb friction)
// We solve them one at a time, each by an impulse that fixes its own constraint given the current
// velocities. Fixing one disturbs its neighbours, so we sweep over all of them several times and the
// velocities converge. The key to making repeated sweeps correct is to clamp the ACCUMULATED impulse
// (sum over all sweeps), not each sweep's increment: a contact may push, but never pull.
#include "body.hpp"
#include "collision.hpp"

#include <vector>

namespace phys {

struct SolverSettings {
    int iterations = 10;
    // Fraction of the remaining overlap removed per second, as a velocity (Baumgarte stabilisation).
    // Too high adds energy and makes resting bodies jitter; too low lets them sink.
    Real baumgarte = static_cast<Real>(0.2);
    // Overlap tolerated without correction. Keeping a sliver of contact stops resting bodies from
    // alternating between touching and separated every frame, which would lose their contacts.
    Real slop = static_cast<Real>(0.005);
    // Closing speeds below this do not bounce (otherwise a body settling on the ground would
    // chatter forever on ever smaller bounces).
    Real restitution_threshold = 1;
};

// One colliding pair found by the narrow phase: indices into the world's body array.
struct ContactPair {
    int a = 0, b = 0;
    Manifold manifold;
};

// Changes the velocities of the bodies in `contacts` so none are moving into each other, bouncing
// and applying friction as their materials dictate. `dt` converts overlap into a correcting speed.
void solve_contacts(std::vector<Body>& bodies, const std::vector<ContactPair>& contacts, Real dt,
                    const SolverSettings& settings);

}  // namespace phys
