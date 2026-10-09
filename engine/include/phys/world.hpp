#pragma once
// The simulation loop. A World owns the bodies and advances them one fixed step at a time.
#include "broadphase.hpp"
#include "solver.hpp"

#include <cstdint>
#include <vector>

namespace phys {

enum class BroadphaseKind {
    BruteForce,     // every pair: the O(n^2) reference
    SweepAndPrune,  // sorted by x, persistent order
    DynamicTree,    // padded AABBs in a balanced tree (the default)
};

// What the last step did, for the demo's readout and for tests.
struct StepStats {
    std::size_t candidate_pairs = 0;     // pairs the broad phase handed to the narrow phase
    std::size_t contacts = 0;            // pairs the narrow phase confirmed touching
    std::uint64_t broadphase_tests = 0;  // box comparisons / tree nodes visited by the broad phase
};

class World {
public:
    Vec2 gravity{0, static_cast<Real>(-9.81)};
    SolverSettings solver;
    BroadphaseKind broadphase = BroadphaseKind::DynamicTree;
    std::vector<Body> bodies;

    // Returns the body's index. Indices stay valid as long as bodies are only ever appended.
    int add(const Body& body) {
        bodies.push_back(body);
        return static_cast<int>(bodies.size()) - 1;
    }

    // Keeps only the first `count` bodies (e.g. to delete everything added after the level was built).
    // Use this rather than resizing `bodies` directly: it also drops the remembered contacts, whose
    // body indices would otherwise point at the wrong bodies and warm start them with stale impulses.
    void truncate(size_t count) {
        bodies.resize(count);
        contacts_.clear();
    }

    // One step:
    //   1. broad phase: AABB overlaps give candidate pairs; filters drop pairs that must not collide
    //   2. narrow phase on the candidates, then match each contact point to last step's (by id) so it
    //      can inherit that step's impulse
    //   3. apply gravity/forces to velocities
    //   4. change velocities so contacts stop approaching (bounce + friction), warm started
    //   5. move bodies using those velocities
    void step(Real dt);

    // Contacts found during the most recent step (for debug drawing and tests).
    const std::vector<ContactPair>& contacts() const { return contacts_; }
    const StepStats& stats() const { return stats_; }
    // The tree broad phase's structure, for drawing. Only meaningful once a step has used it.
    const DynamicTree& broadphase_tree() const { return tree_.tree(); }

private:
    std::vector<ContactPair> contacts_;
    StepStats stats_;
    SweepAndPrune sap_;
    TreeBroadphase tree_;
    std::vector<AABB> boxes_;                // scratch, reused every step
    std::vector<std::uint8_t> immovable_;
};

}  // namespace phys
