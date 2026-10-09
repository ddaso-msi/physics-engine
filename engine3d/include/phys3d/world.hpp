#pragma once
// Stage 10, step 4: the 3D simulation loop.
//
// Deliberately the minimum that makes things collide and stack: every pair is tested (after a cheap
// bounding-sphere rejection), there are no joints, no sleeping and no continuous collision yet. Those come
// from the 2D engine's designs in later steps.
#include "solver.hpp"

#include <vector>

namespace phys3d {

class World {
public:
    Vec3 gravity{0, static_cast<Real>(-9.81), 0};
    SolverSettings solver;
    Gyroscopic gyroscopic = Gyroscopic::Implicit;
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

    // One step: find contacts and match them to last step's, apply gravity and forces to velocities,
    // solve the contacts, move the bodies.
    void step(Real dt);

    // Contacts found during the most recent step.
    const std::vector<ContactPair>& contacts() const { return contacts_; }

private:
    std::vector<ContactPair> contacts_;
};

}  // namespace phys3d
