#pragma once
// Win conditions, as plain questions about the simulation state. No input, no drawing: a level asks
// these after each step, and the tests ask them about hand-built worlds.
#include "level.hpp"

#include <span>

namespace puzzle {

// How many of `bodies` (indices into w.bodies) have their centre of mass below height `y`.
int count_below(const World& w, std::span<const int> bodies, Real y);

// Is the body's centre of mass inside the axis-aligned region [lo, hi]?
bool inside(const World& w, int body, Vec2 lo, Vec2 hi);

// Are all of `bodies` slower than both limits (m/s and rad/s)? A sleeping body counts as slow.
bool all_slow(const World& w, std::span<const int> bodies, Real linear = 0.1f, Real angular = 0.2f);

// The highest point reached by any of `bodies`' outlines. kFloorTop if the list is empty.
Real top_of(const World& w, std::span<const int> bodies);

// "True for long enough": a condition has to hold on every step for `required` simulated seconds.
// One step where it does not hold starts the count again, so a single lucky frame never passes.
struct Hold {
    double held = 0;
    // Feed it the condition after each step. Returns true once it has held for `required` seconds.
    bool update(bool condition, Real dt, double required) {
        held = condition ? held + static_cast<double>(dt) : 0;
        return held >= required;
    }
};

}  // namespace puzzle
