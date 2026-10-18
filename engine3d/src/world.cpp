#include <phys3d/world.hpp>

#include <algorithm>
#include <chrono>
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
std::uint64_t pair_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) | static_cast<std::uint32_t>(b);
}

Quat blend(const Quat& q0, const Quat& q1, Real t) { return (q0 * (1 - t) + q1 * t).normalized(); }

}  // namespace

void World::step(Real dt) {
    const size_t n = bodies.size();
    // Where the time goes, for the benchmark: each lap() adds the time since the last one to a phase.
    using Clock = std::chrono::steady_clock;
    Clock::time_point mark = Clock::now();
    auto lap = [&](double& phase_ms) {
        const Clock::time_point now = Clock::now();
        phase_ms = std::chrono::duration<double, std::milli>(now - mark).count();
        mark = now;
    };

    // 0. Anything that disturbs a sleeper wakes it: a force or torque, or a velocity someone set.
    for (Body& b : bodies) {
        if (b.type != BodyType::Dynamic || b.awake) continue;
        if (!allow_sleep || b.force.length_sq() != 0 || b.torque.length_sq() != 0 || b.vel.length_sq() != 0 || b.w.length_sq() != 0)
            b.wake();
    }
    // A motor told to turn, with torque to do it, must be able to start a sleeping mechanism: commanding
    // it touches no velocity and no force.
    for (const Joint& j : joints) {
        if ((j.type == JointType::Hinge || j.type == JointType::Slider) && j.enable_motor && j.max_motor > 0 && j.motor_speed != 0) {
            if (j.a >= 0) bodies[static_cast<size_t>(j.a)].wake();
            if (j.b >= 0) bodies[static_cast<size_t>(j.b)].wake();
        }
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
        const int rebuilds = tree_.rebuilds();
        candidates = tree_.find_pairs(boxes_, immovable_);
        stats_.broadphase_tests = tree_.last_tests();
        stats_.broadphase_rebuilt = tree_.rebuilds() != rebuilds;
    }
    stats_.candidate_pairs = candidates.size();

    // Bodies held together by a joint do not collide with each other (unless the joint says they may).
    no_collide_.clear();
    for (const Joint& j : joints)
        if (!j.collide_connected && j.a >= 0 && j.b >= 0) no_collide_.insert(pair_key(j.a, j.b));

    lap(stats_.ms_broadphase);

    // 2. Narrow phase. A pair where neither body is being simulated cannot change, so it is not recomputed...
    //    Each pair's GJK starts from where it ended last step, if the pair was a candidate then too. Both
    //    lists are sorted by pair, so finding last step's entry is a walk along it, not a search.
    gjk_cache_next_.clear();
    size_t cached = 0;
    stats_.gjk_pairs = 0;
    stats_.gjk_iterations = 0;
    for (const IndexPair& pair : candidates) {
        const Body& a = bodies[static_cast<size_t>(pair.a)];
        const Body& b = bodies[static_cast<size_t>(pair.b)];
        if (!a.is_active() && !b.is_active()) continue;
        const std::uint64_t key = pair_key(pair.a, pair.b);
        if (!no_collide_.empty() && no_collide_.count(key)) continue;
        while (cached < gjk_cache_.size() && gjk_cache_[cached].key < key) ++cached;
        GjkCache gjk;
        if (cached < gjk_cache_.size() && gjk_cache_[cached].key == key) gjk = gjk_cache_[cached].gjk;
        gjk.iterations = 0;
        ContactPair contact;
        const bool touching = collide_pair(a, b, contact.manifold, &gjk);
        if (gjk.iterations > 0) {  // (zero: this pair never reached GJK, the per-shape code handled it)
            ++stats_.gjk_pairs;
            stats_.gjk_iterations += static_cast<std::uint64_t>(gjk.iterations);
        }
        gjk_cache_next_.push_back({key, gjk});
        if (touching) {
            contact.a = pair.a;
            contact.b = pair.b;
            contacts_.push_back(contact);
        }
    }
    std::swap(gjk_cache_, gjk_cache_next_);
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

    lap(stats_.ms_narrowphase);

    // 3. Islands. A contact or a joint between two dynamic bodies links them; static bodies link nothing (a floor of
    //    separate piles is not one island).
    UnionFind islands(n);
    auto dynamic = [&](int i) { return bodies[static_cast<size_t>(i)].type == BodyType::Dynamic; };
    for (const ContactPair& c : contacts_)
        if (c.manifold.count > 0 && dynamic(c.a) && dynamic(c.b)) islands.unite(c.a, c.b);
    for (const Joint& j : joints)
        if (j.a >= 0 && j.b >= 0 && dynamic(j.a) && dynamic(j.b)) islands.unite(j.a, j.b);

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

    lap(stats_.ms_islands);

    // 4-6. Forces, constraints, movement.
    for (Body& b : bodies)
        if (b.awake) b.integrate_velocity(dt, gravity, gyroscopic);
    solve_constraints(bodies, contacts_, joints, dt, solver, &joint_solver_);
    lap(stats_.ms_solver);
    prev_pos_.resize(n);
    prev_q_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        prev_pos_[i] = bodies[i].pos;
        prev_q_[i] = bodies[i].q;
    }
    for (Body& b : bodies)
        if (b.awake) b.integrate_position(dt);

    lap(stats_.ms_movement);

    // 7. Continuous collision.
    stats_.ccd_swept = 0;
    stats_.ccd_hits = 0;
    if (continuous) continuous_collision();

    lap(stats_.ms_sweep);

    // 7b. The joint position pass, after the sweep (step 30): the sweep pulls a fast body back along its path with no
    // thought for its joints, and done before it the pass left what the sweep opened (8 cm, once, on a ragdoll
    // tumbling down steps) for the next step. (Counted with the movement.) But the pass knows nothing of walls either:
    // closing a joint can pull a body the sweep had just stopped into a thin wall, from where it went on through. So
    // the sweep is run again, from where each body began the step to where the joints have left it: what the pass
    // put through a wall is stopped short of it (and a joint held open against a wall stays so, for the contacts to
    // settle). A body already stopped the first time is not moved again: it is just short of the wall.
    solve_joint_positions(bodies, joints, solver, fixed_contacts());
    double joints_ms = 0;
    lap(joints_ms);
    stats_.ms_movement += joints_ms;
    if (continuous && !joints.empty()) {
        const size_t swept = stats_.ccd_swept;
        continuous_collision();
        stats_.ccd_swept = swept;  // (the same bodies, counted once)
        // What the second sweep stopped, it stopped short of a wall the pass had pulled it into, and the joints
        // it pulled on are open again. (A ragdoll's arm stopped by a step: the pass closed its shoulder mostly by
        // moving the arm, the lighter part, into the step, and the sweep put the arm back. Up to 4 cm open, a few
        // times a run, and much of the time again the next step.) So the pass is run once more with those bodies
        // held where the sweep left them, as if fixed: the rest of what is joined to them comes to them instead.
        // Only for the joints of their islands; the rest have had their pass. Then swept again, for what that
        // moved. (Two balls fired at a wall are both stopped by the first sweep, the reason not to hold what that
        // one stops: the pair's joint would stay as open as the sweep left it.)
        if (std::find(stopped_.begin(), stopped_.end(), 1) != stopped_.end()) {
            held_.clear();
            island_held_.assign(n, 0);
            for (size_t i = 0; i < n; ++i)
                if (stopped_[i]) {
                    Body& b = bodies[i];
                    held_.push_back({i, b.inv_mass, b.inv_inertia_body});
                    b.type = BodyType::Static;
                    b.inv_mass = 0;
                    b.inv_inertia_body = Mat3::diagonal(0, 0, 0);
                    island_held_[static_cast<size_t>(islands.find(static_cast<int>(i)))] = 1;
                }
            auto in_held_island = [&](int body) {
                return body >= 0 && island_held_[static_cast<size_t>(islands.find(body))] != 0;
            };
            held_joints_.clear();
            for (const Joint& j : joints)  // (copies: the pass only moves bodies)
                if (in_held_island(j.a) || in_held_island(j.b)) held_joints_.push_back(j);
            solve_joint_positions(bodies, held_joints_, solver, fixed_contacts());
            for (const Held& h : held_) {
                Body& b = bodies[h.index];
                b.type = BodyType::Dynamic;
                b.inv_mass = h.inv_mass;
                b.inv_inertia_body = h.inv_inertia_body;
            }
            for (size_t i = 0; i < n; ++i) island_held_[i] = in_held_island(static_cast<int>(i));  // (now per body)
            continuous_collision(&island_held_);
            stats_.ccd_swept = swept;
        }
        double again_ms = 0;
        lap(again_ms);
        stats_.ms_sweep += again_ms;
    }

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
    lap(stats_.ms_sleep);
}

// This step's contacts of jointed bodies with fixed ones, for the joint position pass (step 33): rows along the normal
// out of the fixed body, each with how far the pass may move the body along it: no deeper than the slop, or than it
// is already. Each pair the narrow phase found touching is looked at again where the bodies are now, after their move
// and the sweep. Still overlapping, a row per point of its contact. Apart (the velocity solve will often have lifted a
// body off what it landed on), one row between the closest points, which the pass may close by the gap and the slop:
// with no row there, a ball dragged along a ledge was pulled 6 cm into it by the pass, and the sweep let it stay 2.5 cm
// deeper every step until it came out underneath. Apart by more than kNear, none. (Rows from the narrow phase's
// contact instead held a ball that had slid off the end of a ledge to the ledge's top as if it went on, and its joint
// was left open, by 9 cm at 60 m/s.)
std::span<FixedContact> World::fixed_contacts() {
    fixed_contacts_.clear();
    if (joints.empty()) return fixed_contacts_;
    jointed_.assign(bodies.size(), 0);
    for (const Joint& j : joints)
        for (const int b : {j.a, j.b})
            if (b >= 0) jointed_[static_cast<size_t>(b)] = 1;
    for (const ContactPair& pair : contacts_) {
        if (pair.dormant) continue;
        const Body& a = bodies[static_cast<size_t>(pair.a)];
        const Body& b = bodies[static_cast<size_t>(pair.b)];
        const bool a_fixed = a.type != BodyType::Dynamic, b_fixed = b.type != BodyType::Dynamic;
        if (a_fixed == b_fixed) continue;
        const int moving = a_fixed ? pair.b : pair.a;
        const Body& body = bodies[static_cast<size_t>(moving)];
        if (!jointed_[static_cast<size_t>(moving)] || !body.is_active()) continue;
        constexpr Real kNear = static_cast<Real>(0.25);  // (more than the pass moves a body in one round)
        const ClosestResult apart = closest(Convex::of(a), Convex::of(b), /*with_points=*/true);
        if (apart.distance > kNear) continue;
        if (apart.distance > 0) {
            const Vec3 out = a_fixed ? apart.normal : -apart.normal, at = a_fixed ? apart.point_b : apart.point_a;
            fixed_contacts_.push_back({moving, out, cross(at - body.pos, out), -(apart.distance + solver.slop), nullptr});
            continue;
        }
        Manifold now;
        if (!collide_pair(a, b, now)) continue;
        const Vec3 out = a_fixed ? now.normal : -now.normal;
        for (int k = 0; k < now.count; ++k) {
            const ContactPoint& cp = now.points[k];
            fixed_contacts_.push_back({moving, out, cross(cp.point - body.pos, out), std::min(Real(0), cp.depth - solver.slop), nullptr});
        }
    }
    return fixed_contacts_;
}

// The 2D engine's sweep (see phys/world.cpp for the full reasoning), in 3D. A body that moves or turns far
// enough in one step to jump over a thin obstacle is tested at poses spaced no further apart than its own
// inscribed radius, which a convex body at least that wide cannot step over.
//   - Starting clear of the wall: the first overlapping sample brackets the time of impact, bisection
//     refines it, and the body is put at the overlapping end, just barely inside, so next step's solver sees
//     a contact and responds.
//   - Starting already in contact: its centre may advance no more than a quarter of the inscribed radius
//     further into the wall in one step, measured along the direction it entered by. If it is held there,
//     the part of its velocity that was carrying it in is reversed like any bounce.
//
// That second rule is the third version. Capping the overlap depth collide() reports does not work: for a
// plate thinner than the wall the depth stops growing as the plate slices in. Following the body's farthest
// point along the entry direction does not work either: a thin capsule that arrives end-on and turns side-on
// as it goes in reaches no deeper with its farthest point (the end swings back as fast as the body advances)
// while all of it passes through. And holding the pose without touching the velocity is not enough: when
// the one contact point the solver has is swinging AWAY from the wall (a long body spinning hard), the solver
// does nothing, and the body crept in a quarter radius per step until it was through.
void World::continuous_collision(const std::vector<std::uint8_t>* only) {
    const size_t n = bodies.size();
    stopped_.assign(n, 0);
    for (size_t i = 0; i < n; ++i) {
        Body& body = bodies[i];
        if (!body.is_active() || (only && !(*only)[i])) continue;
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
            *touching = i < wall_index ? collide_pair(probe, wall, m) : collide_pair(wall, probe, m);
            return *touching ? m.depth : Real(0);
        };

        const int steps = std::clamp(static_cast<int>(std::ceil(motion / radius)), 1, 512);
        Real earliest = 2;  // > 1 means "no hit"
        bool held = false;  // the earliest stop was a body already in contact reaching its cap...
        Vec3 held_into;     // ...moving this way into the wall
        Real held_bounce = 0;
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
                touching = i < s ? collide_pair(probe, wall, start) : collide_pair(wall, probe, start);
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
                if (hi < earliest) earliest = hi, held = false;
            } else {
                // Already touching. d points from the body into the wall; the centre travels in a straight
                // line, so the fraction of the step at which it has gone a further depth_cap that way is a
                // simple ratio.
                const Vec3 d = i < s ? start.normal : -start.normal;
                const Real travel = dot(p1 - p0, d);
                if (travel <= depth_cap) continue;
                const Real stop = depth_cap / travel;
                if (stop < earliest) {
                    earliest = stop;
                    held = true;
                    held_into = d;
                    held_bounce = std::max(body.restitution, wall.restitution);
                }
            }
        }
        if (earliest <= 1) {
            body.pos = p0 + (p1 - p0) * earliest;
            body.q = blend(q0, q1, earliest);
            ++stats_.ccd_hits;
            stopped_[i] = 1;
            if (held) {
                // Held at its cap: it is being driven into the wall and the contact solver did nothing
                // about it this step, or it would not be here. Treat it as the impact it is.
                const Real closing = dot(body.vel, held_into);
                if (closing > 0) body.vel -= held_into * ((1 + (closing > solver.restitution_threshold ? held_bounce : Real(0))) * closing);
            }
        }
    }
}

}  // namespace phys3d
