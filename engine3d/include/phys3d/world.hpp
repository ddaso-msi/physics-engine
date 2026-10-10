#pragma once
// Stage 10, step 4: the 3D simulation loop.
//
// Contacts, joints, sleeping and continuous collision.
#include "broadphase.hpp"
#include "solver.hpp"

#include <cstdint>
#include <unordered_set>
#include <vector>

namespace phys3d {

enum class BroadphaseKind {
    BruteForce,   // every pair of AABBs: the O(n^2) reference
    DynamicTree,  // padded AABBs in a balanced tree (the default)
};

// What the last step did, for the demo's readout and for tests.
struct StepStats {
    std::size_t candidate_pairs = 0;     // pairs the broad phase handed to the narrow phase
    std::size_t contacts = 0;            // contacts that took part in solving (dormant ones excluded)
    std::uint64_t broadphase_tests = 0;  // box comparisons / tree nodes visited
    std::size_t awake_bodies = 0;        // dynamic bodies being simulated (not asleep)
    std::size_t islands = 0;             // groups of dynamic bodies connected by contacts or joints
    std::size_t ccd_swept = 0;           // bodies moving fast enough to be swept by the continuous pass
    std::size_t ccd_hits = 0;            // of those, how many it stopped at a static body
};

class World {
public:
    Vec3 gravity{0, static_cast<Real>(-9.81), 0};
    SolverSettings solver;
    Gyroscopic gyroscopic = Gyroscopic::Implicit;
    BroadphaseKind broadphase = BroadphaseKind::DynamicTree;
    std::vector<Body> bodies;
    std::vector<Joint> joints;

    // Sleeping, as in the 2D engine. An island (bodies linked by contacts and joints) goes to sleep once every body in
    // it has been slower than both thresholds for `time_to_sleep` seconds. A sleeping body is not
    // integrated or solved; anything that touches it, or sets its velocity or a force on it, wakes the
    // whole island.
    bool allow_sleep = true;
    Real time_to_sleep = static_cast<Real>(0.5);
    Real sleep_linear_speed = static_cast<Real>(0.01);    // m/s
    Real sleep_angular_speed = static_cast<Real>(0.035);  // rad/s, about 2 degrees per second

    // Continuous collision against static bodies, as in the 2D engine. A body that moves or turns more
    // than half its inscribed radius in a step is swept along its path; if it would have gone into or
    // through a static body it is placed where it first touches, and the contact solver deals with it
    // next step. (Fast dynamic-vs-dynamic tunnelling is not handled.)
    bool continuous = true;

    // Returns the body's index.
    int add(const Body& body) {
        bodies.push_back(body);
        return static_cast<int>(bodies.size()) - 1;
    }
    // Joints connect bodies by index, and are solved together with the contacts every step.
    int add_joint(const Joint& joint) {
        joints.push_back(joint);
        return static_cast<int>(joints.size()) - 1;
    }
    // Wakes a body (and, on the next step, everything touching or jointed to it). Call this after moving a body by hand;
    // changing its velocity or applying a force wakes it automatically.
    void wake(int index) { bodies[static_cast<size_t>(index)].wake(); }

    // Keeps only the first `count` bodies, and forgets the contacts (their indices would be stale) and the
    // joints that referred to a removed body.
    void truncate(size_t count) {
        bodies.resize(count);
        contacts_.clear();
        const int limit = static_cast<int>(count);
        std::erase_if(joints, [&](const Joint& j) { return j.a >= limit || j.b >= limit; });
    }

    // One step:
    //   0. wake sleepers that were disturbed (a force, a velocity, a motor told to turn)
    //   1. broad phase: AABB overlaps give candidate pairs
    //   2. narrow phase on the candidates (skipping jointed pairs and pairs where both bodies sleep), then match each contact
    //      point to last step's by id so it can inherit that step's impulse
    //   3. build islands from contacts and joints; an awake body touching a sleeping island wakes it
    //   4. apply gravity and forces to velocities
    //   5. solve the joints and contacts
    //   6. move the bodies, then correct joint position error
    //   7. continuous collision: stop fast bodies that went into or through a static body
    //   8. put islands that have been still long enough to sleep
    void step(Real dt);

    // Contacts found during the most recent step. Includes dormant ones.
    const std::vector<ContactPair>& contacts() const { return contacts_; }
    const StepStats& stats() const { return stats_; }
    const DynamicTree& broadphase_tree() const { return tree_.tree(); }

private:
    void continuous_collision();

    std::vector<ContactPair> contacts_;
    std::vector<Vec3> prev_pos_;  // poses before this step's movement, for the continuous pass
    std::vector<Quat> prev_q_;
    std::vector<std::uint8_t> island_awake_;
    std::vector<Real> island_still_;
    StepStats stats_;
    TreeBroadphase tree_;
    std::vector<AABB> boxes_;  // scratch, reused every step
    std::vector<std::uint8_t> immovable_;
    std::unordered_set<std::uint64_t> no_collide_;  // body pairs joined by a joint that forbids their collision
};

}  // namespace phys3d
