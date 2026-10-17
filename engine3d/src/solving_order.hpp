#pragma once
// The order constraints are solved in (step 19). Internal to the engine.
//
// Gauss-Seidel takes constraints one at a time, each from the velocities the last one left, and two in a
// row that share a body cannot overlap: the second waits for the first's result. In the order they come
// they nearly always share one: contacts sorted by body pair ((5, 7), (5, 9), (5, 12), ...), the joints of a
// chain link to link. So they are coloured, greedily: each takes the lowest colour that no earlier
// constraint on either of its moving bodies has, and they are solved colour by colour. Within a colour no two
// share a moving body, and the processor can work on several at once. A static body (or index -1, the
// world) does not count, because nothing is ever written to it.
#include <phys3d/body.hpp>

#include <cstdint>
#include <vector>

namespace phys3d {

// The number of colours. The last is shared by everything that finds the others taken (a body in more than 63
// constraints), so within it two constraints may share a body.
constexpr int kColours = 64;

// `bodies_of(i, a, b)` sets the two body indices of constraint i and returns false to leave it out. If `starts` is
// given it receives, for each colour, where its constraints begin in the order, and last where the order ends.
template <class BodiesOf>
std::vector<size_t> solving_order(const std::vector<Body>& bodies, size_t count, BodiesOf bodies_of, std::vector<size_t>* starts = nullptr) {
    std::vector<std::uint64_t> taken(bodies.size(), 0);  // per body: the colours its constraints have so far
    std::vector<int> colour(count, -1);
    int per_colour[kColours] = {};
    auto moves = [&](int i) { return i >= 0 && bodies[static_cast<size_t>(i)].type == BodyType::Dynamic; };
    for (size_t i = 0; i < count; ++i) {
        int a = -1, b = -1;
        if (!bodies_of(i, a, b)) continue;
        const std::uint64_t used = (moves(a) ? taken[static_cast<size_t>(a)] : 0) | (moves(b) ? taken[static_cast<size_t>(b)] : 0);
        int k = 0;
        while (k < kColours - 1 && (used >> k & 1)) ++k;  // (a body in 64 constraints shares the last colour)
        colour[i] = k;
        ++per_colour[k];
        if (moves(a)) taken[static_cast<size_t>(a)] |= std::uint64_t(1) << k;
        if (moves(b)) taken[static_cast<size_t>(b)] |= std::uint64_t(1) << k;
    }
    int start[kColours];
    int total = 0;
    for (int k = 0; k < kColours; ++k) start[k] = total, total += per_colour[k];
    if (starts) {
        starts->assign(start, start + kColours);
        starts->push_back(static_cast<size_t>(total));
    }
    std::vector<size_t> order(static_cast<size_t>(total));
    for (size_t i = 0; i < count; ++i)
        if (colour[i] >= 0) order[static_cast<size_t>(start[colour[i]]++)] = i;
    return order;
}

}  // namespace phys3d
