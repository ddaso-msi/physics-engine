#pragma once
// Win conditions, as plain questions about the simulation state. No input, no drawing: a level asks
// these after each step, and the tests ask them about hand-built worlds.
#include "level.hpp"

#include <span>

namespace puzzle {

// How many of `bodies` (indices into w.bodies) have their centre of mass below height `y`.
int count_below(const World& w, std::span<const int> bodies, Real y);

}  // namespace puzzle
