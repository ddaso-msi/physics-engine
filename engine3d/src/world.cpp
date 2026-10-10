#include <phys3d/world.hpp>

#include <utility>

namespace phys3d {

void World::step(Real dt) {
    const std::vector<ContactPair> previous = std::move(contacts_);
    contacts_.clear();

    // Broad phase. Both algorithms return the same pairs (the tree, a superset) in the same order.
    const size_t n = bodies.size();
    boxes_.resize(n);
    immovable_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        boxes_[i] = compute_aabb(bodies[i]);
        immovable_[i] = bodies[i].inv_mass == 0;
    }
    PairList candidates;
    if (broadphase == BroadphaseKind::BruteForce) {
        candidates = brute_force_pairs(boxes_, immovable_, &stats_.broadphase_tests);
    } else {
        candidates = tree_.find_pairs(boxes_, immovable_);
        stats_.broadphase_tests = tree_.last_tests();
    }
    stats_.candidate_pairs = candidates.size();

    // Narrow phase.
    for (const IndexPair& pair : candidates) {
        ContactPair contact;
        if (collide(bodies[static_cast<size_t>(pair.a)], bodies[static_cast<size_t>(pair.b)], contact.manifold)) {
            contact.a = pair.a;
            contact.b = pair.b;
            contacts_.push_back(contact);
        }
    }
    stats_.contacts = contacts_.size();
    transfer_impulses(previous, contacts_);

    for (Body& b : bodies) b.integrate_velocity(dt, gravity, gyroscopic);
    solve_contacts(bodies, contacts_, dt, solver);
    for (Body& b : bodies) b.integrate_position(dt);
}

}  // namespace phys3d
