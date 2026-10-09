#pragma once
// Stage 7: broad phase. Narrow-phase collision (collision.hpp) is exact but costs real work per pair,
// and testing every pair is O(n^2). The broad phase cheaply discards pairs that cannot touch by
// comparing their AABBs, and hands the narrow phase only the survivors ("candidate pairs").
//
// All three implementations take the same input (one box per body, indexed like World::bodies, plus a
// flag marking immovable bodies) and return pairs sorted by (a, b) with a < b. Sorting makes the
// output order independent of the algorithm, so the simulation is bit-for-bit identical whichever
// broad phase runs. Pairs of two immovable bodies are never reported: they can never respond.
#include "aabb.hpp"
#include "dynamic_tree.hpp"

#include <compare>
#include <cstdint>
#include <span>
#include <vector>

namespace phys {

struct IndexPair {
    int a = 0, b = 0;  // a < b
    auto operator<=>(const IndexPair&) const = default;
};
using PairList = std::vector<IndexPair>;

// Reference implementation: test every pair. O(n^2); every other broad phase must agree with it.
// (`tests`, if given, receives the number of box comparisons made.)
PairList brute_force_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable,
                           std::uint64_t* tests = nullptr);

// Sweep and prune. Keep the bodies sorted by the left edge of their box along x; a body can only
// overlap the bodies after it whose left edge is no further right than its own right edge, so the
// inner loop stops early. The sorted order is kept between calls: bodies barely move per step, so
// re-sorting with insertion sort is close to O(n). Cost is about O(n sqrt n) for evenly spread bodies
// but degrades when many boxes share an x range (a tall stack, a wide floor).
class SweepAndPrune {
public:
    PairList find_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable);
    // Box comparisons in the last call, and shifts the last re-sort needed (for tests / benchmarks).
    std::uint64_t last_tests() const { return tests_; }
    std::uint64_t last_shifts() const { return shifts_; }

private:
    std::vector<int> order_;  // body indices sorted by box.lo.x
    std::uint64_t tests_ = 0, shifts_ = 0;
};

// Dynamic AABB tree broad phase. Each body's box is stored in the tree padded by `margin`; the tree is
// only touched when a body's real box leaves its padded one, so most steps cost a containment test
// per body. Each moving body then queries the tree, about O(log n) each.
class TreeBroadphase {
public:
    explicit TreeBroadphase(Real margin = static_cast<Real>(0.1)) : margin_(margin) {}

    // Returns a superset of brute_force_pairs: bodies whose real boxes overlap are always reported,
    // plus some whose padded boxes overlap. The narrow phase discards the extras.
    PairList find_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable);

    const DynamicTree& tree() const { return tree_; }
    // Number of proxies that had to be moved in the last call (bodies that escaped their padding).
    int last_moves() const { return moves_; }
    std::uint64_t last_tests() const { return tests_; }

private:
    DynamicTree tree_;
    std::vector<int> proxy_;  // body index -> tree proxy id
    Real margin_;
    int moves_ = 0;
    std::uint64_t tests_ = 0;
};

}  // namespace phys
