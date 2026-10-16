#pragma once
// Stage 10, step 6: broad phase for 3D. As in 2D, it compares cheap AABBs and hands the narrow phase only
// the pairs that might touch. Both implementations take one box per body (indexed like World::bodies) and
// a flag marking immovable bodies, and return pairs sorted by (a, b) with a < b, never a pair of two
// immovable bodies. Sorting makes the order independent of the algorithm, so the simulation is bit-for-bit
// the same whichever runs.
#include "aabb.hpp"
#include "dynamic_tree.hpp"

#include <compare>
#include <cstdint>
#include <span>
#include <vector>

namespace phys3d {

struct IndexPair {
    int a = 0, b = 0;  // a < b
    auto operator<=>(const IndexPair&) const = default;
};
using PairList = std::vector<IndexPair>;

// Reference implementation: test every pair. O(n^2). (`tests`, if given, receives the comparison count.)
PairList brute_force_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable, std::uint64_t* tests = nullptr);

// Dynamic AABB tree. Each body's box is stored padded by `margin`; the tree is only touched when a body's
// real box leaves its padded one. The overlapping pairs are then found by one pass of the tree over itself.
// Returns exactly the pairs brute_force_pairs does, in the same order. The tree is laid out afresh when a
// quarter of its leaves have been re-inserted, or when searching it has become a tenth dearer than it was
// just after the last time.
class TreeBroadphase {
public:
    explicit TreeBroadphase(Real margin = static_cast<Real>(0.1)) : margin_(margin) {}

    PairList find_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable);

    const DynamicTree& tree() const { return tree_; }
    int last_moves() const { return moves_; }           // proxies re-inserted in the last call
    int rebuilds() const { return rebuilds_; }          // times the whole tree has been built afresh
    std::uint64_t last_tests() const { return tests_; }  // tree nodes visited in the last call

private:
    DynamicTree tree_;
    std::vector<int> proxy_;  // body index -> tree proxy id
    Real margin_;
    int moves_ = 0;
    std::size_t churn_ = 0;  // leaves created or re-inserted since the tree was last built from scratch
    int rebuilds_ = 0;
    std::uint64_t tests_ = 0;
    double fresh_cost_ = 0;  // nodes visited per pair and leaf by the first search after the last rebuild
    bool worn_ = false;      // the last search cost enough more than that to call for a rebuild
};

}  // namespace phys3d
