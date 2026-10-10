#include <phys3d/world.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

namespace phys3d {

namespace {

// Islands are the connected components of "touching", which is exactly what union-find computes.
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

// Orientation a fraction t of the way from q0 to q1. q1 was reached by integrating forward from q0, so the
// two are on the same side of the quaternion sphere and a straight blend follows the way the body actually
// turned (taking "the short way round" instead would be wrong for a turn of more than half a revolution).
// A normalised blend is not perfectly even in angle, which does not matter: it only has to pass through
// every pose in between.
Quat blend(const Quat& q0, const Quat& q1, Real t) { return (q0 * (1 - t) + q1 * t).normalized(); }

// The point of the body's shape farthest along unit direction d, if the body were at (pos, q).
Vec3 support_point(const Body& b, Vec3 pos, const Quat& q, Vec3 d) {
    if (b.shape.type == Shape::Type::Sphere) return pos + d * b.shape.radius;
    const Mat3 r = to_mat3(q);
    const Vec3 axis[3] = {r.cx, r.cy, r.cz};
    Vec3 p = pos;
    for (int k = 0; k < 3; ++k) p += axis[k] * (dot(axis[k], d) < 0 ? -b.shape.half_extents[k] : b.shape.half_extents[k]);
    return p;
}

}  // namespace

void World::step(Real dt) {
    const size_t n = bodies.size();

    // 0. Anything that disturbs a sleeper wakes it: a force or torque, or a velocity someone set.
    for (Body& b : bodies) {
        if (b.type != BodyType::Dynamic || b.awake) continue;
        if (!allow_sleep || b.force.length_sq() != 0 || b.torque.length_sq() != 0 || b.vel.length_sq() != 0 || b.w.length_sq() != 0)
            b.wake();
    }

    std::vector<ContactPair> previous = std::move(contacts_);
    contacts_.clear();

    // 1. Broad phase. Both algorithms return the same pairs (the tree, a superset) in the same order.
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

    // 2. Narrow phase. A pair where neither body is being simulated cannot change, so it is not recomputed...
    for (const IndexPair& pair : candidates) {
        const Body& a = bodies[static_cast<size_t>(pair.a)];
        const Body& b = bodies[static_cast<size_t>(pair.b)];
        if (!a.is_active() && !b.is_active()) continue;
        ContactPair contact;
        if (collide(a, b, contact.manifold)) {
            contact.a = pair.a;
            contact.b = pair.b;
            contacts_.push_back(contact);
        }
    }
    // ...but a contact between two sleeping bodies (or a sleeper and the ground) is remembered, with its
    // impulses, so that waking the pile warm starts it instead of letting it settle from scratch.
    for (ContactPair& c : previous) {
        if (static_cast<size_t>(c.a) >= n || static_cast<size_t>(c.b) >= n) continue;
        if (!bodies[static_cast<size_t>(c.a)].is_active() && !bodies[static_cast<size_t>(c.b)].is_active()) {
            c.dormant = true;
            contacts_.push_back(c);
        }
    }
    transfer_impulses(previous, contacts_);

    // 3. Islands. A contact between two dynamic bodies links them; static bodies link nothing (a floor of
    //    separate piles is not one island).
    UnionFind islands(n);
    auto dynamic = [&](int i) { return bodies[static_cast<size_t>(i)].type == BodyType::Dynamic; };
    for (const ContactPair& c : contacts_)
        if (c.manifold.count > 0 && dynamic(c.a) && dynamic(c.b)) islands.unite(c.a, c.b);

    island_awake_.assign(n, 0);
    stats_.islands = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!dynamic(static_cast<int>(i))) continue;
        const size_t root = static_cast<size_t>(islands.find(static_cast<int>(i)));
        if (root == i) ++stats_.islands;
        if (bodies[i].awake) island_awake_[root] = 1;
    }
    // A sleeper in an island that has an awake body is part of the action now: wake it.
    for (size_t i = 0; i < n; ++i)
        if (dynamic(static_cast<int>(i)) && !bodies[i].awake && island_awake_[static_cast<size_t>(islands.find(static_cast<int>(i)))])
            bodies[i].wake();
    for (ContactPair& c : contacts_)
        if (c.dormant && (bodies[static_cast<size_t>(c.a)].is_active() || bodies[static_cast<size_t>(c.b)].is_active()))
            c.dormant = false;
    stats_.contacts = static_cast<size_t>(std::count_if(contacts_.begin(), contacts_.end(), [](const ContactPair& c) { return !c.dormant; }));

    // 4-6. Forces, contacts, movement.
    for (Body& b : bodies)
        if (b.awake) b.integrate_velocity(dt, gravity, gyroscopic);
    solve_contacts(bodies, contacts_, dt, solver);
    prev_pos_.resize(n);
    prev_q_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        prev_pos_[i] = bodies[i].pos;
        prev_q_[i] = bodies[i].q;
    }
    for (Body& b : bodies)
        if (b.awake) b.integrate_position(dt);

    // 7. Continuous collision.
    stats_.ccd_swept = 0;
    stats_.ccd_hits = 0;
    if (continuous) continuous_collision();

    // 8. Sleep. A body is "still" while it is slower than both thresholds; an island sleeps when its most
    //    restless member has been still long enough.
    if (allow_sleep) {
        const Real lin2 = sleep_linear_speed * sleep_linear_speed, ang2 = sleep_angular_speed * sleep_angular_speed;
        for (Body& b : bodies) {
            if (!b.is_active()) continue;
            const bool still = b.vel.length_sq() <= lin2 && b.w.length_sq() <= ang2;
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
                bodies[i].w = {};
                bodies[i].sleep_time = 0;
            }
        }
    }
    stats_.awake_bodies = static_cast<size_t>(std::count_if(bodies.begin(), bodies.end(), [](const Body& b) { return b.is_active(); }));
}

// The 2D engine's sweep (see phys/world.cpp for the full reasoning), in 3D. A body that moves or turns far
// enough in one step to jump over a thin obstacle is tested at poses spaced no further apart than its own
// inscribed radius, which a convex body at least that wide cannot step over.
//   - Starting clear of the wall: the first overlapping sample brackets the time of impact, bisection
//     refines it, and the body is put at the overlapping end, just barely inside, so next step's solver sees
//     a contact and responds.
//   - Starting already in contact: its deepest point may advance no more than a quarter of the inscribed
//     radius further into the wall in one step. "Deeper" is measured along the direction it entered by,
//     as the travel of its farthest point that way. (Capping the overlap depth collide() reports does not
//     work: for a plate thinner than the wall that depth stops growing as the plate slices in, and starts
//     measuring the way out the far side.)
void World::continuous_collision() {
    const size_t n = bodies.size();
    for (size_t i = 0; i < n; ++i) {
        Body& body = bodies[i];
        if (!body.is_active()) continue;
        const Vec3 p0 = prev_pos_[i], p1 = body.pos;
        const Quat q0 = prev_q_[i], q1 = body.q;
        const Real radius = body.inscribed_radius(), reach = body.bounding_radius();
        const Real turned = rotation_angle(q1 * conjugate(q0));
        const Real motion = distance(p0, p1) + turned * reach;  // how far its outline travels
        if (radius <= 0 || motion <= static_cast<Real>(0.5) * radius) continue;  // slow enough: ordinary contacts suffice
        ++stats_.ccd_swept;

        const Vec3 pad{reach, reach, reach};
        const AABB corridor{Vec3{std::min(p0.x, p1.x), std::min(p0.y, p1.y), std::min(p0.z, p1.z)} - pad,
                            Vec3{std::max(p0.x, p1.x), std::max(p0.y, p1.y), std::max(p0.z, p1.z)} + pad};
        const Real depth_cap = static_cast<Real>(0.25) * radius;

        // How deep the body is inside static body `wall_index` at fraction t of its path (0 if clear). It
        // asks exactly the narrow phase's question, lower index first, so the two can never disagree about
        // a razor-thin overlap.
        auto depth_at = [&](Real t, size_t wall_index, bool* touching) {
            const Body& wall = bodies[wall_index];
            Body probe = body;
            probe.pos = p0 + (p1 - p0) * t;
            probe.q = blend(q0, q1, t);
            Manifold m;
            *touching = i < wall_index ? collide(probe, wall, m) : collide(wall, probe, m);
            return *touching ? m.depth : Real(0);
        };

        const int steps = std::clamp(static_cast<int>(std::ceil(motion / radius)), 1, 512);
        Real earliest = 2;  // > 1 means "no hit"
        for (size_t s = 0; s < n; ++s) {
            const Body& wall = bodies[s];
            if (s == i || wall.inv_mass != 0) continue;  // only static bodies
            if (!overlap(corridor, boxes_[s])) continue;  // boxes_ still holds this step's AABBs; statics do not move

            bool touching;
            Manifold start;
            {
                Body probe = body;
                probe.pos = p0;
                probe.q = q0;
                touching = i < s ? collide(probe, wall, start) : collide(wall, probe, start);
            }
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
            } else {
                // Already touching. d points from the body into the wall; follow the body's farthest point
                // along d and stop it once that point has gone a further depth_cap in.
                const Vec3 d = i < s ? start.normal : -start.normal;
                const Real start_reach = dot(support_point(body, p0, q0, d), d);
                auto advance = [&](Real t) {
                    return dot(support_point(body, p0 + (p1 - p0) * t, blend(q0, q1, t), d), d) - start_reach;
                };
                int first = -1;
                for (int k = 1; k <= steps && first < 0; ++k)
                    if (advance(static_cast<Real>(k) / static_cast<Real>(steps)) > depth_cap) first = k;
                if (first < 0) continue;
                Real lo = static_cast<Real>(first - 1) / static_cast<Real>(steps), hi = static_cast<Real>(first) / static_cast<Real>(steps);
                for (int it = 0; it < 16; ++it) {
                    const Real mid = (lo + hi) * static_cast<Real>(0.5);
                    (advance(mid) > depth_cap ? hi : lo) = mid;
                }
                earliest = std::min(earliest, lo);  // the last pose still within the cap
            }
        }
        if (earliest <= 1) {
            body.pos = p0 + (p1 - p0) * earliest;
            body.q = blend(q0, q1, earliest);
            ++stats_.ccd_hits;
        }
    }
}

}  // namespace phys3d
