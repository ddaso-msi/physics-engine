#include <phys/world.hpp>

#include <algorithm>
#include <utility>

namespace phys {

namespace {
// Order-independent key for a pair of body indices.
std::uint64_t pair_key(int a, int b) {
    const auto lo = static_cast<std::uint64_t>(static_cast<std::uint32_t>(std::min(a, b)));
    const auto hi = static_cast<std::uint64_t>(static_cast<std::uint32_t>(std::max(a, b)));
    return (lo << 32) | hi;
}
}  // namespace

void World::step(Real dt) {
    const std::vector<ContactPair> previous = std::move(contacts_);
    contacts_.clear();

    // Broad phase. Everything here is a function of the body list, so all three algorithms return the
    // same pairs (or, for the tree, a superset) in the same order.
    const size_t n = bodies.size();
    boxes_.resize(n);
    immovable_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        boxes_[i] = compute_aabb(bodies[i]);
        immovable_[i] = bodies[i].inv_mass == 0;
    }
    PairList candidates;
    switch (broadphase) {
        case BroadphaseKind::BruteForce:
            candidates = brute_force_pairs(boxes_, immovable_, &stats_.broadphase_tests);
            break;
        case BroadphaseKind::SweepAndPrune:
            candidates = sap_.find_pairs(boxes_, immovable_);
            stats_.broadphase_tests = sap_.last_tests();
            break;
        case BroadphaseKind::DynamicTree:
            candidates = tree_.find_pairs(boxes_, immovable_);
            stats_.broadphase_tests = tree_.last_tests();
            break;
    }
    stats_.candidate_pairs = candidates.size();

    // Bodies held together by a joint do not collide with each other (unless the joint says they may).
    no_collide_.clear();
    for (const Joint& j : joints)
        if (!j.collide_connected && j.a >= 0 && j.b >= 0) no_collide_.insert(pair_key(j.a, j.b));

    // Narrow phase.
    for (const IndexPair& pair : candidates) {
        const Body& a = bodies[static_cast<size_t>(pair.a)];
        const Body& b = bodies[static_cast<size_t>(pair.b)];
        if (!should_collide(a, b)) continue;
        if (!no_collide_.empty() && no_collide_.count(pair_key(pair.a, pair.b))) continue;
        ContactPair contact;
        if (collide(a, b, contact.manifold)) {
            contact.a = pair.a;
            contact.b = pair.b;
            contacts_.push_back(contact);
        }
    }
    stats_.contacts = contacts_.size();
    transfer_impulses(previous, contacts_);

    for (Body& b : bodies) b.integrate_velocity(dt, gravity);
    solve_constraints(bodies, contacts_, joints, dt, solver);
    for (Body& b : bodies) b.integrate_position(dt);
    solve_joint_positions(bodies, joints, solver);
}

}  // namespace phys
