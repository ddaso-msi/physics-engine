#pragma once
// The simulation loop. A World owns the bodies and advances them one fixed step at a time.
#include "broadphase.hpp"
#include "solver.hpp"

#include <cstdint>
#include <unordered_set>
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
    std::vector<Joint> joints;

    // Returns the body's index. Indices stay valid as long as bodies are only ever appended.
    int add(const Body& body) {
        bodies.push_back(body);
        return static_cast<int>(bodies.size()) - 1;
    }

    // Joints connect bodies by index, and are solved together with the contacts every step.
    int add_joint(const Joint& joint) {
        joints.push_back(joint);
        return static_cast<int>(joints.size()) - 1;
    }
    // Removes one joint. Later joints shift down by one index.
    void remove_joint(int index) { joints.erase(joints.begin() + index); }

    // Keeps only the first `count` bodies (e.g. to delete everything added after the level was built).
    // Use this rather than resizing `bodies` directly: it also drops the remembered contacts, whose
    // body indices would otherwise point at the wrong bodies and warm start them with stale impulses,
    // and the joints that referred to a removed body.
    void truncate(size_t count) {
        bodies.resize(count);
        contacts_.clear();
        const int limit = static_cast<int>(count);
        std::erase_if(joints, [&](const Joint& j) { return j.a >= limit || j.b >= limit; });
    }

    // One step:
    //   1. broad phase: AABB overlaps give candidate pairs; filters drop pairs that must not collide
    //   2. narrow phase on the candidates (skipping jointed pairs), then match each contact point to
    //      last step's (by id) so it can inherit that step's impulse
    //   3. apply gravity/forces to velocities
    //   4. change velocities so joints hold and contacts stop approaching (bounce + friction), warm started
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
    std::unordered_set<std::uint64_t> no_collide_;  // body pairs joined by a joint that forbids their collision
};

}  // namespace phys
