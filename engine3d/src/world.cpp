#include <phys3d/world.hpp>

#include <utility>

namespace phys3d {

namespace {
// Radius of a sphere around the centre of mass that contains the whole shape.
Real bounding_radius(const Body& b) {
    return b.shape.type == Shape::Type::Sphere ? b.shape.radius : b.shape.half_extents.length();
}
}  // namespace

void World::step(Real dt) {
    const std::vector<ContactPair> previous = std::move(contacts_);
    contacts_.clear();

    // Every pair: O(n^2), with a bounding-sphere test to skip the far-apart ones cheaply.
    const int n = static_cast<int>(bodies.size());
    for (int i = 0; i < n; ++i) {
        const Body& a = bodies[static_cast<size_t>(i)];
        const Real ra = bounding_radius(a);
        for (int j = i + 1; j < n; ++j) {
            const Body& b = bodies[static_cast<size_t>(j)];
            if (a.inv_mass == 0 && b.inv_mass == 0) continue;  // two immovable bodies never interact
            const Real reach = ra + bounding_radius(b);
            if ((b.pos - a.pos).length_sq() > reach * reach) continue;
            ContactPair pair;
            if (collide(a, b, pair.manifold)) {
                pair.a = i;
                pair.b = j;
                contacts_.push_back(pair);
            }
        }
    }
    transfer_impulses(previous, contacts_);

    for (Body& b : bodies) b.integrate_velocity(dt, gravity, gyroscopic);
    solve_contacts(bodies, contacts_, dt, solver);
    for (Body& b : bodies) b.integrate_position(dt);
}

}  // namespace phys3d
