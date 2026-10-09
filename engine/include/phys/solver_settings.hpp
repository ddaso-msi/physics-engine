#pragma once
// Tuning knobs shared by the contact solver and the joint solver.
#include "math.hpp"

namespace phys {

struct SolverSettings {
    int iterations = 10;
    // Fraction of the remaining overlap removed per second, as a velocity (Baumgarte stabilisation),
    // for CONTACTS. Too high adds energy and makes resting bodies jitter; too low lets them sink.
    // Joints do not use it: they correct position error in a separate pass (joint_position_iterations).
    Real baumgarte = static_cast<Real>(0.2);
    // Passes of the joint position correction that runs after bodies have moved.
    int joint_position_iterations = 4;
    // Overlap tolerated without correction. Keeping a sliver of contact stops resting bodies from
    // alternating between touching and separated every frame, which would lose their contacts.
    Real slop = static_cast<Real>(0.005);
    // Closing speeds below this do not bounce (otherwise a body settling on the ground would
    // chatter forever on ever smaller bounces).
    Real restitution_threshold = 1;
    // Start each step from last step's impulses (found by World via the contact point ids) instead
    // of from zero. A resting contact needs nearly the same impulse every frame, so this begins the
    // sweeps almost at the answer, and it lets a stack's weight propagate down over a few frames
    // instead of needing many sweeps in a single one. Joints keep their impulses between steps too.
    bool warm_starting = true;
};

}  // namespace phys
