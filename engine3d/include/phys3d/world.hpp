#pragma once
// Stage 10, step 4: the 3D simulation loop.
//
// Contacts only so far: there are no joints, no sleeping and no continuous collision yet. Those come from
// the 2D engine's designs in later steps.
#include "broadphase.hpp"
#include "solver.hpp"

#include <cstdint>

#include <vector>

namespace phys3d {

enum class BroadphaseKind {
    BruteForce,   // every pair of AABBs: the O(n^2) reference
    DynamicTree,  // padded AABBs in a balanced tree (the default)
};

// What the last step did, for the demo's readout and for tests.
struct StepStats {
    std::size_t candidate_pairs = 0;     // pairs the broad phase handed to the narrow phase
    std::size_t contacts = 0;            // pairs the narrow phase confirmed touching
    std::uint64_t broadphase_tests = 0;  // box comparisons / tree nodes visited
};

class World {
public:
    Vec3 gravity{0, static_cast<Real>(-9.81), 0};
    SolverSettings solver;
    Gyroscopic gyroscopic = Gyroscopic::Implicit;
    BroadphaseKind broadphase = BroadphaseKind::DynamicTree;
    std::vector<Body> bodies;

    // Returns the body's index.
    int add(const Body& body) {
        bodies.push_back(body);
        return static_cast<int>(bodies.size()) - 1;
    }
    // Keeps only the first `count` bodies, and forgets the contacts (their indices would be stale).
    void truncate(size_t count) {
        bodies.resize(count);
        contacts_.clear();
    }

    // One step: broad phase (AABB overlaps give candidate pairs), narrow phase on the candidates, match
    // contacts to last step's, apply gravity and forces to velocities, solve the contacts, move the bodies.
    void step(Real dt);

    // Contacts found during the most recent step.
    const std::vector<ContactPair>& contacts() const { return contacts_; }
    const StepStats& stats() const { return stats_; }
    const DynamicTree& broadphase_tree() const { return tree_.tree(); }

private:
    std::vector<ContactPair> contacts_;
    StepStats stats_;
    TreeBroadphase tree_;
    std::vector<AABB> boxes_;  // scratch, reused every step
    std::vector<std::uint8_t> immovable_;
};

}  // namespace phys3d
