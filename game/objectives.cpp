#include "objectives.hpp"

namespace puzzle {

int count_below(const World& w, std::span<const int> bodies, Real y) {
    int n = 0;
    for (int i : bodies)
        if (w.bodies[static_cast<size_t>(i)].pos.y < y) ++n;
    return n;
}

}  // namespace puzzle
