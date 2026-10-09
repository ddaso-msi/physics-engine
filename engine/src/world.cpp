#include <phys/world.hpp>

namespace phys {

void World::step(Real dt) {
    contacts_.clear();
    const int n = static_cast<int>(bodies.size());
    // Every pair: O(n^2). Stage 7 replaces this with a broad phase.
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const Body& a = bodies[static_cast<size_t>(i)];
            const Body& b = bodies[static_cast<size_t>(j)];
            if (a.inv_mass == 0 && b.inv_mass == 0) continue;  // two immovable bodies never interact
            ContactPair pair;
            if (collide(a, b, pair.manifold)) {
                pair.a = i;
                pair.b = j;
                contacts_.push_back(pair);
            }
        }
    }

    for (Body& b : bodies) b.integrate_velocity(dt, gravity);
    solve_contacts(bodies, contacts_, dt, solver);
    for (Body& b : bodies) b.integrate_position(dt);
}

}  // namespace phys
