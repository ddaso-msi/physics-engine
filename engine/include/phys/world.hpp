#pragma once
// Stage 5: the simulation loop. A World owns the bodies and advances them one fixed step at a time.
#include "solver.hpp"

#include <vector>

namespace phys {

class World {
public:
    Vec2 gravity{0, static_cast<Real>(-9.81)};
    SolverSettings solver;
    std::vector<Body> bodies;

    // Returns the body's index. Indices stay valid as long as bodies are only ever appended.
    int add(const Body& body) {
        bodies.push_back(body);
        return static_cast<int>(bodies.size()) - 1;
    }

    // One step:
    //   1. find contacts at the current positions, and match each point to last step's (by id) so
    //      it can inherit that step's impulse
    //   2. apply gravity/forces to velocities
    //   3. change velocities so contacts stop approaching (bounce + friction), warm started
    //   4. move bodies using those velocities
    void step(Real dt);

    // Keeps only the first `count` bodies (e.g. to delete everything added after the level was built).
    // Use this rather than resizing `bodies` directly: it also drops the remembered contacts, whose
    // body indices would otherwise point at the wrong bodies and warm start them with stale impulses.
    void truncate(size_t count) {
        bodies.resize(count);
        contacts_.clear();
    }

    // Contacts found during the most recent step (for debug drawing and tests).
    const std::vector<ContactPair>& contacts() const { return contacts_; }

private:
    std::vector<ContactPair> contacts_;
};

}  // namespace phys
