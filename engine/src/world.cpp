#include <phys/world.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

namespace phys {

namespace {

// Order-independent key for a pair of body indices.
std::uint64_t pair_key(int a, int b) {
    const auto lo = static_cast<std::uint64_t>(static_cast<std::uint32_t>(std::min(a, b)));
    const auto hi = static_cast<std::uint64_t>(static_cast<std::uint32_t>(std::max(a, b)));
    return (lo << 32) | hi;
}

// Islands are the connected components of "touching or jointed", which is exactly what union-find computes.
struct UnionFind {
    std::vector<int> parent;
    explicit UnionFind(size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    int find(int x) {
        while (parent[static_cast<size_t>(x)] != x) {
            parent[static_cast<size_t>(x)] = parent[static_cast<size_t>(parent[static_cast<size_t>(x)])];  // path halving
            x = parent[static_cast<size_t>(x)];
        }
        return x;
    }
    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[static_cast<size_t>(std::max(a, b))] = std::min(a, b);
    }
};

// Farthest a body's shape reaches from its centre of mass.
Real bounding_radius(const Body& b) {
    if (b.shape.type == Shape::Type::Circle) return b.shape.circle.radius;
    Real r2 = 0;
    for (int i = 0; i < b.shape.polygon.count; ++i) r2 = std::max(r2, b.shape.polygon.vertices[i].length_sq());
    return std::sqrt(r2);
}

}  // namespace

void World::step(Real dt) {
    const size_t n = bodies.size();

    // 0. Anything that disturbs a sleeper wakes it: a force or torque, a velocity someone set, a mouse
    //    joint holding it (the user is dragging it, so it must stay live).
    for (Body& b : bodies) {
        if (b.type != BodyType::Dynamic || b.awake) continue;
        if (!allow_sleep || b.force.length_sq() != 0 || b.torque != 0 || b.vel.length_sq() != 0 || b.w != 0) b.wake();
    }
    for (const Joint& j : joints) {
        // A motor told to turn, with torque to do it, has to be able to start a sleeping mechanism: nothing
        // about commanding it touches a velocity or a force, so without this a sleeping car ignores its throttle.
        const bool driven = (j.type == JointType::Revolute || j.type == JointType::Prismatic) && j.enable_motor &&
                            j.max_motor > 0 && j.motor_speed != 0;
        if (j.type == JointType::Mouse || driven) {
            if (j.a >= 0) bodies[static_cast<size_t>(j.a)].wake();
            if (j.b >= 0) bodies[static_cast<size_t>(j.b)].wake();
        }
    }

    std::vector<ContactPair> previous = std::move(contacts_);
    contacts_.clear();

    // 1. Broad phase. Everything here is a function of the body list, so all three algorithms return the
    //    same pairs (or, for the tree, a superset) in the same order.
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

    // 2. Narrow phase. A pair where neither body is being simulated cannot change, so it is not recomputed.
    for (const IndexPair& pair : candidates) {
        const Body& a = bodies[static_cast<size_t>(pair.a)];
        const Body& b = bodies[static_cast<size_t>(pair.b)];
        if (!a.is_active() && !b.is_active()) continue;
        if (!should_collide(a, b)) continue;
        if (!no_collide_.empty() && no_collide_.count(pair_key(pair.a, pair.b))) continue;
        ContactPair contact;
        if (collide(a, b, contact.manifold)) {
            contact.a = pair.a;
            contact.b = pair.b;
            contacts_.push_back(contact);
        }
    }
    // ...but a contact between two sleeping bodies (or a sleeper and the ground) must be remembered, with
    // its impulses, so that waking the pile warm starts it instead of letting it settle from scratch.
    for (ContactPair& c : previous) {
        if (static_cast<size_t>(c.a) >= n || static_cast<size_t>(c.b) >= n) continue;
        if (!bodies[static_cast<size_t>(c.a)].is_active() && !bodies[static_cast<size_t>(c.b)].is_active()) {
            c.dormant = true;
            contacts_.push_back(c);
        }
    }
    transfer_impulses(previous, contacts_);

    // 3. Islands. Contacts and joints between two dynamic bodies link them; static bodies do not (a whole
    //    floor of resting piles is not one island).
    UnionFind islands(n);
    auto dynamic = [&](int i) { return bodies[static_cast<size_t>(i)].type == BodyType::Dynamic; };
    for (const ContactPair& c : contacts_)
        if (c.manifold.count > 0 && dynamic(c.a) && dynamic(c.b)) islands.unite(c.a, c.b);
    for (const Joint& j : joints)
        if (j.a >= 0 && j.b >= 0 && dynamic(j.a) && dynamic(j.b)) islands.unite(j.a, j.b);

    island_awake_.assign(n, 0);
    island_asleep_.assign(n, 0);
    stats_.islands = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!dynamic(static_cast<int>(i))) continue;
        const size_t root = static_cast<size_t>(islands.find(static_cast<int>(i)));
        if (root == i) ++stats_.islands;
        (bodies[i].awake ? island_awake_ : island_asleep_)[root] = 1;
    }
    // A sleeper in an island that has an awake body is part of the action now: wake it.
    for (size_t i = 0; i < n; ++i)
        if (dynamic(static_cast<int>(i)) && !bodies[i].awake && island_awake_[static_cast<size_t>(islands.find(static_cast<int>(i)))])
            bodies[i].wake();
    for (ContactPair& c : contacts_)
        if (c.dormant && (bodies[static_cast<size_t>(c.a)].is_active() || bodies[static_cast<size_t>(c.b)].is_active()))
            c.dormant = false;
    stats_.contacts = static_cast<size_t>(std::count_if(contacts_.begin(), contacts_.end(), [](const ContactPair& c) { return !c.dormant; }));

    // 4-6. Forces, constraints, movement.
    for (Body& b : bodies)
        if (b.awake) b.integrate_velocity(dt, gravity);
    solve_constraints(bodies, contacts_, joints, dt, solver);
    prev_pos_.resize(n);
    prev_angle_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        prev_pos_[i] = bodies[i].pos;
        prev_angle_[i] = bodies[i].angle;
    }
    for (Body& b : bodies)
        if (b.awake) b.integrate_position(dt);
    solve_joint_positions(bodies, joints, solver);

    // 7. Continuous collision.
    stats_.ccd_swept = 0;
    stats_.ccd_hits = 0;
    if (continuous) continuous_collision(dt);

    // 8. Sleep. A body is "still" while it is slower than both thresholds; an island sleeps when its
    //    most restless member has been still long enough.
    if (allow_sleep) {
        const Real lin2 = sleep_linear_speed * sleep_linear_speed;
        for (Body& b : bodies) {
            if (!b.is_active()) continue;
            const bool still = b.vel.length_sq() <= lin2 && std::fabs(b.w) <= sleep_angular_speed;
            b.sleep_time = (still && b.allow_sleep) ? b.sleep_time + dt : 0;
        }
        island_still_.assign(n, std::numeric_limits<Real>::infinity());
        for (size_t i = 0; i < n; ++i)
            if (bodies[i].is_active()) {
                const size_t root = static_cast<size_t>(islands.find(static_cast<int>(i)));
                island_still_[root] = std::min(island_still_[root], bodies[i].sleep_time);
            }
        for (size_t i = 0; i < n; ++i) {
            if (!bodies[i].is_active()) continue;
            if (island_still_[static_cast<size_t>(islands.find(static_cast<int>(i)))] >= time_to_sleep) {
                bodies[i].awake = false;
                bodies[i].vel = {};
                bodies[i].w = 0;
                bodies[i].sleep_time = 0;
            }
        }
    }
    stats_.awake_bodies = static_cast<size_t>(
        std::count_if(bodies.begin(), bodies.end(), [](const Body& b) { return b.is_active(); }));
}

// A body that moves farther in one step than it is wide can jump clean over a thin wall without ever
// overlapping it at either end of the step, so the narrow phase never sees the wall. Sweep such bodies
// along their path instead. "Moves" includes turning: a fast spinner sweeps its corners through space too.
//
// Test the body at poses spaced no further apart than its own inscribed radius. A convex body that is at
// least 2r wide and takes steps of at most r must overlap any obstacle it crosses at one of the samples,
// however thin the obstacle is.
//   - Starting clear of the wall: the first overlapping sample brackets the time of impact, bisection
//     refines it, and the body is put at the overlapping end of the bracket: just barely inside the wall,
//     so next step the ordinary contact solver sees a contact and responds (bounce, friction).
//   - Starting already in contact (just after an impact, say): the solver cannot always hold a body that
//     is still moving fast, because an off-centre hit can turn most of the momentum into a swing around the
//     contact point. So it may sink no deeper than a quarter of its inscribed radius in one step; the
//     sweep stops it where it reaches that depth. A quarter radius is shallower than the way out on
//     either side of any wall, so the contact still pushes it back out the side it came in from.
void World::continuous_collision(Real dt) {
    const size_t n = bodies.size();
    for (size_t i = 0; i < n; ++i) {
        Body& body = bodies[i];
        if (!body.is_active()) continue;
        const Vec2 p0 = prev_pos_[i], p1 = body.pos;
        // The turn it actually made this step. NOT the difference of the two stored angles: those are kept
        // wrapped into (-pi, pi], so beyond half a revolution per step their difference is the short way
        // round, the wrong way, and the sweep would test poses the body never passed through.
        const Real a0 = prev_angle_[i], da = body.w * dt;
        const Real radius = body.inscribed_radius();
        const Real reach = bounding_radius(body);
        const Real motion = distance(p0, p1) + std::fabs(da) * reach;  // how far its outline travels
        if (radius <= 0 || motion <= static_cast<Real>(0.5) * radius) continue;  // slow enough: ordinary contacts suffice
        ++stats_.ccd_swept;

        // Everything the body can touch on the way lies within its bounding radius of the path.
        const AABB corridor{{std::min(p0.x, p1.x) - reach, std::min(p0.y, p1.y) - reach},
                            {std::max(p0.x, p1.x) + reach, std::max(p0.y, p1.y) + reach}};
        const Real depth_cap = static_cast<Real>(0.25) * radius;

        // How deep the body is inside static body `wall_index` at fraction t of its path (0 if clear). Asks
        // exactly the question the narrow phase will ask next step, including which body is passed first
        // (the lower index). At the razor-thin overlaps the search homes in on, rounding can make
        // collide(a, b) and collide(b, a) disagree; if the sweep said "touching" but the narrow phase said
        // "not touching" the body would go on through and never be stopped.
        auto depth_at = [&](Real t, size_t wall_index, bool* touching) {
            const Body& wall = bodies[wall_index];
            Body probe = body;
            probe.pos = p0 + (p1 - p0) * t;
            probe.set_angle(a0 + da * t);
            Manifold m;
            *touching = i < wall_index ? collide(probe, wall, m) : collide(wall, probe, m);
            return *touching ? m.depth : Real(0);
        };

        const int steps = std::clamp(static_cast<int>(std::ceil(motion / radius)), 1, 512);
        Real earliest = 2;  // > 1 means "no hit"
        for (size_t s = 0; s < n; ++s) {
            const Body& wall = bodies[s];
            if (s == i || wall.inv_mass != 0) continue;  // only static bodies
            if (!should_collide(body, wall) || !overlap(corridor, compute_aabb(wall))) continue;

            bool touching;
            const Real depth0 = depth_at(0, s, &touching);
            if (!touching) {
                // Clear at the start: find where it first overlaps.
                int first = -1;
                for (int k = 1; k <= steps && first < 0; ++k) {
                    depth_at(static_cast<Real>(k) / static_cast<Real>(steps), s, &touching);
                    if (touching) first = k;
                }
                if (first < 0) continue;
                Real lo = static_cast<Real>(first - 1) / static_cast<Real>(steps), hi = static_cast<Real>(first) / static_cast<Real>(steps);
                for (int it = 0; it < 16; ++it) {
                    const Real mid = (lo + hi) * static_cast<Real>(0.5);
                    depth_at(mid, s, &touching);
                    (touching ? hi : lo) = mid;
                }
                earliest = std::min(earliest, hi);
            } else if (depth0 <= depth_cap) {
                // Already touching and not too deep: find where it would sink past the cap.
                int first = -1;
                for (int k = 1; k <= steps && first < 0; ++k)
                    if (depth_at(static_cast<Real>(k) / static_cast<Real>(steps), s, &touching) > depth_cap) first = k;
                if (first < 0) continue;
                Real lo = static_cast<Real>(first - 1) / static_cast<Real>(steps), hi = static_cast<Real>(first) / static_cast<Real>(steps);
                for (int it = 0; it < 16; ++it) {
                    const Real mid = (lo + hi) * static_cast<Real>(0.5);
                    (depth_at(mid, s, &touching) > depth_cap ? hi : lo) = mid;
                }
                earliest = std::min(earliest, lo);  // the last pose still within the cap
            }
            // else: already deeper than the cap at the start; the solver has to deal with it.
        }
        if (earliest <= 1) {
            body.pos = p0 + (p1 - p0) * earliest;
            body.set_angle(a0 + da * earliest);
            ++stats_.ccd_hits;
        }
    }
}

}  // namespace phys
