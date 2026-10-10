#include "test.hpp"
#include <phys3d/world.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

using namespace phys3d;

namespace {

constexpr double kG = 9.81;
constexpr Real kDt = 1.0f / 120.0f;

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
};

Body floor_box() {
    Body f = Body::fixed_box({30, 0.5f, 30}, {0, -0.5f, 0});
    f.restitution = 0;
    return f;
}
Body crate(Vec3 p, Real half = 0.5f) {
    Body b = Body::solid_box({half, half, half}, 1, p);
    b.restitution = 0;
    return b;
}
Body ball(Real r, Vec3 p) {
    Body b = Body::solid_sphere(r, 1, p);
    b.restitution = 0;
    return b;
}
World with_floor() {
    World w;
    w.add(floor_box());
    return w;
}
void run(World& w, int steps) {
    for (int i = 0; i < steps; ++i) w.step(kDt);
}
void add_pyramid(World& w, int base) {  // base x base, then (base-1) x (base-1), ...
    for (int layer = 0; layer < base; ++layer) {
        const int m = base - layer;
        for (int i = 0; i < m; ++i)
            for (int k = 0; k < m; ++k)
                w.add(crate({(static_cast<Real>(i) - static_cast<Real>(m - 1) * 0.5f) * 1.02f, 0.5f + static_cast<Real>(layer) * 1.002f,
                             (static_cast<Real>(k) - static_cast<Real>(m - 1) * 0.5f) * 1.02f}));
    }
}
bool same_pose(const Body& a, const Body& b) {
    return a.pos.x == b.pos.x && a.pos.y == b.pos.y && a.pos.z == b.pos.z && a.q.x == b.q.x && a.q.y == b.q.y && a.q.z == b.q.z && a.q.w == b.q.w;
}
// A thin static wall across x = 8: 4 cm thick, far taller and wider than any shot can miss.
void add_thin_wall(World& w) { w.add(Body::fixed_box({0.02f, 30, 30}, {8, 5, 0})); }

}  // namespace

// ---- sleeping ------------------------------------------------------------------------------------

TEST(a_resting_crate_falls_asleep_and_then_does_not_move_at_all_in_3d) {
    World w = with_floor();
    const size_t c = static_cast<size_t>(w.add(crate({0, 1.5f, 0})));
    run(w, 360);
    CHECK(!w.bodies[c].awake);
    CHECK(w.stats().awake_bodies == 0);
    CHECK(w.bodies[c].vel.length_sq() == 0 && w.bodies[c].w.length_sq() == 0);

    const Body before = w.bodies[c];
    run(w, 240);
    CHECK(same_pose(w.bodies[c], before));  // bit for bit

    // Its contact with the floor is remembered (dormant) along with the impulses that held it up.
    CHECK(!w.contacts().empty());
    double total = 0;
    for (const ContactPair& cp : w.contacts()) {
        CHECK(cp.dormant);
        for (int k = 0; k < cp.manifold.count; ++k) total += double(cp.manifold.points[k].normal_impulse);
    }
    CHECK_NEAR(total, kG * double(kDt), 0.05 * kG * double(kDt));
    CHECK(w.stats().contacts == 0);  // nothing was solved
}

TEST(a_tumbling_body_does_not_sleep_and_slow_drift_is_told_from_rest) {
    World tumbling;
    tumbling.gravity = {};
    Body brick = Body::solid_box({1.0f, 0.6f, 0.2f}, 1, {0, 0, 0});
    brick.w = {3, 2, 1};
    tumbling.add(brick);
    run(tumbling, 480);
    CHECK(tumbling.bodies[0].awake);

    auto asleep_after_1s = [](Vec3 vel, Vec3 spin) {
        World w;
        w.gravity = {};
        Body b = ball(0.2f, {0, 0, 0});
        b.vel = vel;
        b.w = spin;
        w.add(b);
        run(w, 120);
        return !w.bodies[0].awake;
    };
    CHECK(!asleep_after_1s({0.02f, 0, 0}, {}));  // above 0.01 m/s: drifting
    CHECK(asleep_after_1s({0.005f, 0, 0}, {}));  // below it: asleep after half a second
    CHECK(!asleep_after_1s({}, {0, 0.06f, 0}));  // turning faster than 0.035 rad/s: not at rest
    CHECK(asleep_after_1s({}, {0, 0.02f, 0}));
}

TEST(an_island_sleeps_and_wakes_as_one_in_3d) {
    World w = with_floor();
    for (int i = 0; i < 4; ++i) w.add(crate({0, 0.5f + static_cast<Real>(i) * 1.001f, 0}));
    int steps = 0;
    for (; steps < 20 * 120; ++steps) {
        w.step(kDt);
        int asleep = 0;
        for (size_t i = 1; i <= 4; ++i) asleep += !w.bodies[i].awake;
        CHECK(asleep == 0 || asleep == 4);  // all at once, never crate by crate
        if (asleep == 4) break;
    }
    CHECK(steps < 20 * 120);

    const size_t dropped = static_cast<size_t>(w.add(ball(0.3f, {0, 7, 0})));
    bool woke = false;
    for (int i = 0; i < 240 && !woke; ++i) {
        w.step(kDt);
        woke = w.bodies[1].awake || w.bodies[2].awake || w.bodies[3].awake || w.bodies[4].awake;
    }
    CHECK(woke);
    for (size_t i = 1; i <= 4; ++i) CHECK(w.bodies[i].awake);  // the whole tower in the same step
    CHECK(w.bodies[dropped].awake);
}

TEST(a_sleeping_pyramid_wakes_warm_started_when_poked) {
    World w = with_floor();
    add_pyramid(w, 3);  // 9 + 4 + 1 crates
    run(w, 6 * 120);
    CHECK(w.stats().awake_bodies == 0);
    std::vector<Vec3> before;
    for (const Body& b : w.bodies) before.push_back(b.pos);
    const size_t top = w.bodies.size() - 1;
    w.bodies[top].apply_impulse_at({0.05f, 0, 0}, w.bodies[top].pos);
    w.step(kDt);
    CHECK(w.bodies[top].awake);
    run(w, 120);
    double worst_below = 0;
    for (size_t i = 1; i < top; ++i) worst_below = std::max(worst_below, double(distance(w.bodies[i].pos, before[i])));
    CHECK(worst_below < 0.02);  // no jump, no sag
}

TEST(velocity_force_or_torque_wake_a_sleeper_and_a_teleport_needs_wake) {
    World w = with_floor();
    const size_t a = static_cast<size_t>(w.add(crate({-4, 0.5f, 0})));
    const size_t b = static_cast<size_t>(w.add(crate({0, 0.5f, 0})));
    const size_t c = static_cast<size_t>(w.add(crate({4, 0.5f, 0})));
    const size_t d = static_cast<size_t>(w.add(crate({8, 0.5f, 0})));
    run(w, 240);
    CHECK(w.stats().awake_bodies == 0);

    w.bodies[a].vel = {0, 0, 2};
    w.bodies[b].apply_force({200, 0, 0});
    w.bodies[c].apply_torque({0, 40, 0});
    w.bodies[d].pos.y = 6;  // moved by hand, nobody told the world
    run(w, 30);
    CHECK(w.bodies[a].pos.z > 0.1f);
    CHECK(w.bodies[b].pos.x > 0.1f);
    CHECK(std::fabs(rotate(w.bodies[c].q, {1, 0, 0}).z) > 0.01f);  // it has turned about the vertical
    CHECK(w.bodies[d].pos.y == 6);                                // still frozen in mid air
    w.wake(static_cast<int>(d));
    run(w, 30);
    CHECK(w.bodies[d].pos.y < 5.9f);
}

TEST(islands_are_counted_and_the_floor_does_not_join_them_in_3d) {
    World apart = with_floor();
    apart.add(crate({-3, 0.5f, 0}));
    apart.add(crate({3, 0.5f, 0}));
    apart.step(kDt);
    CHECK(apart.stats().islands == 2);

    World touching = with_floor();
    touching.add(crate({0, 0.5f, 0}));
    touching.add(crate({0, 0.5f, 0.99f}));
    touching.step(kDt);
    CHECK(touching.stats().islands == 1);
}

TEST(a_body_that_may_not_sleep_keeps_its_island_awake_in_3d) {
    World w = with_floor();
    w.add(crate({0, 0.5f, 0}));
    Body top = crate({0, 1.501f, 0});
    top.allow_sleep = false;
    w.add(top);
    run(w, 600);
    CHECK(w.bodies[1].awake && w.bodies[2].awake);

    World off = with_floor();
    off.allow_sleep = false;
    off.add(crate({0, 0.5f, 0}));
    run(off, 600);
    CHECK(off.bodies[1].awake);
}

TEST(sleeping_saves_work_without_changing_the_result_in_3d) {
    World sleepy = with_floor(), awake = with_floor();
    awake.allow_sleep = false;
    add_pyramid(sleepy, 3);
    add_pyramid(awake, 3);
    run(sleepy, 8 * 120);
    run(awake, 8 * 120);
    CHECK(sleepy.stats().awake_bodies == 0 && awake.stats().awake_bodies == 14);
    CHECK(sleepy.stats().contacts == 0 && awake.stats().contacts > 14);
    double worst = 0;
    for (size_t i = 0; i < sleepy.bodies.size(); ++i) worst = std::max(worst, double(distance(sleepy.bodies[i].pos, awake.bodies[i].pos)));
    CHECK(worst < 0.02);

    sleepy.truncate(1);  // and truncating a world full of sleepers is safe
    CHECK(sleepy.contacts().empty());
    sleepy.add(crate({0, 3, 0}));
    run(sleepy, 240);
    CHECK_NEAR(sleepy.bodies[1].pos.y, 0.5, 0.02);
}

// ---- continuous collision ------------------------------------------------------------------------

TEST(a_fast_ball_stops_at_a_thin_wall_only_with_continuous_collision_in_3d) {
    auto shot = [](bool continuous) {
        World w;
        w.gravity = {};
        w.continuous = continuous;
        add_thin_wall(w);
        Body bullet = ball(0.15f, {2, 5, 0});
        bullet.restitution = 0.5f;
        bullet.vel = {60, 0, 0};  // 0.5 m per step: more than the wall and ball together
        w.add(bullet);
        run(w, 120);
        return w.bodies[1];
    };
    CHECK(shot(false).pos.x > 8.5f);  // sailed straight through
    const Body stopped = shot(true);
    CHECK(stopped.pos.x < 8.0f);
    CHECK_NEAR(stopped.vel.x, -30, 3);  // bounced back at e x the arrival speed
}

TEST(continuous_collision_works_for_any_shot_in_3d) {
    int through_with = 0, through_without = 0, hits = 0;
    for (int trial = 0; trial < 150; ++trial) {
        Lcg rng;
        rng.s = static_cast<std::uint32_t>(trial) * 7919u + 3;
        const Real speed = rng.range(40, 120);
        const Vec3 vel{speed, rng.range(-25, 25), rng.range(-25, 25)};
        const Vec3 spin{rng.range(-12, 12), rng.range(-12, 12), rng.range(-12, 12)};
        const Vec3 start{3, rng.range(2, 8), rng.range(-3, 3)};
        const Quat turn = Quat::from_axis_angle({rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)}, rng.range(0, 3));
        for (int mode = 0; mode < 2; ++mode) {
            World w;
            w.gravity = {};
            w.continuous = mode == 1;
            add_thin_wall(w);
            Body bullet = trial % 2 == 0 ? Body::solid_box({0.12f, 0.2f, 0.16f}, 1, start, turn) : Body::solid_sphere(0.12f, 1, start);
            bullet.vel = vel;
            bullet.w = spin;
            w.add(bullet);
            for (int s = 0; s < 60; ++s) {
                w.step(kDt);
                if (mode == 1) hits += static_cast<int>(w.stats().ccd_hits);
            }
            (mode == 1 ? through_with : through_without) += w.bodies[1].pos.x > 8.0f;
        }
    }
    CHECK(through_with == 0);     // nothing ever gets through
    CHECK(through_without > 90);  // and without it most shots do
    CHECK(hits >= 150);
}

TEST(a_fast_fall_onto_a_thin_floor_does_not_tunnel_in_3d) {
    World w;
    Body plate = Body::fixed_box({20, 0.02f, 20}, {0, 0, 0});
    plate.restitution = 0;
    w.add(plate);
    Body b = ball(0.2f, {0, 6, 0});
    b.vel = {0, -150, 0};
    w.add(b);
    run(w, 120);
    CHECK(w.bodies[1].pos.y > 0.1f);
}

// Turning counts as motion. A long bar spinning in place moves its ends 33 cm a step though its centre
// does not move at all, so it must be swept; the same bar turning slowly, or a settled pile, must not be.
TEST(only_fast_movers_and_fast_spinners_are_swept) {
    auto swept = [](Vec3 spin) {
        World w;
        w.gravity = {};
        Body bar = Body::solid_box({1.0f, 0.05f, 0.05f}, 1, {0, 0, 0});
        bar.w = spin;
        w.add(bar);
        w.step(kDt);
        return w.stats().ccd_swept;
    };
    CHECK(swept({0, 40, 0}) == 1);
    CHECK(swept({0, 0.5f, 0}) == 0);

    World pile = with_floor();
    pile.allow_sleep = false;  // keep every crate live, so "nothing swept" cannot just mean "everything asleep"
    add_pyramid(pile, 3);
    run(pile, 240);
    CHECK(pile.stats().awake_bodies == 14);
    CHECK(pile.stats().ccd_swept == 0);
}

// A bar spinning fast about the vertical beside a thin post: between two steps its tip passes through
// where the post is, though it overlaps the post at neither step. With the sweep it is stopped at the post.
TEST(a_fast_spinner_is_stopped_at_a_post) {
    auto furthest_turn = [](bool continuous) {
        World w;
        w.gravity = {};
        w.continuous = continuous;
        const Real post_angle = 0.5f, reach = 0.9f;  // the post stands on the circle the bar's tip sweeps
        w.add(Body::fixed_box({0.04f, 0.5f, 0.04f}, {reach * std::cos(post_angle), 0, -reach * std::sin(post_angle)}));
        Body bar = Body::solid_box({1.0f, 0.05f, 0.05f}, 1, {0, 0, 0});
        bar.w = {0, 40, 0};  // about +y: its +x end swings toward -z
        bar.restitution = 0;
        w.add(bar);
        w.step(kDt);  // 0.33 rad: not yet at the post
        w.step(kDt);  // 0.67 rad: past the post, unless something stopped it
        const Vec3 tip = rotate(w.bodies[1].q, {1, 0, 0});
        return std::pair<Real, std::size_t>(std::atan2(-tip.z, tip.x), w.stats().ccd_hits);
    };
    const auto without = furthest_turn(false);
    CHECK(without.first > 0.6f);  // swung straight through
    const auto with = furthest_turn(true);
    CHECK(with.first < 0.5f);     // held at the post
    CHECK(with.second == 1);
}

TEST(continuous_collision_leaves_slow_bodies_alone_in_3d) {
    World a = with_floor(), b = with_floor();
    b.continuous = false;
    add_pyramid(a, 3);
    add_pyramid(b, 3);
    for (int i = 0; i < 3; ++i) {
        a.add(ball(0.3f, {-1 + static_cast<Real>(i), 6, 0.2f}));
        b.add(ball(0.3f, {-1 + static_cast<Real>(i), 6, 0.2f}));
    }
    std::size_t hits = 0;
    for (int s = 0; s < 240; ++s) {
        a.step(kDt);
        b.step(kDt);
        hits += a.stats().ccd_hits;
    }
    CHECK(hits == 0);
    for (size_t i = 0; i < a.bodies.size(); ++i) CHECK(same_pose(a.bodies[i], b.bodies[i]));  // bit for bit
}

// ---- regressions found by a 3000-shot stress run ---------------------------------------------------

namespace {

// The stress run's shot number `trial`: a sphere or a box of random size, speed (30 to 200 m/s), spin (up
// to 30 rad/s per axis) and orientation, fired at a wall 1 to 10 cm thick.
World stress_shot(int trial) {
    Lcg rng;
    rng.s = static_cast<std::uint32_t>(trial) * 104729u + 11;
    const Real speed = rng.range(30, 200), thick = rng.range(0.01f, 0.1f), size = rng.range(0.08f, 0.4f);
    const Vec3 vel{speed, rng.range(-40, 40), rng.range(-40, 40)};
    const Vec3 spin{rng.range(-30, 30), rng.range(-30, 30), rng.range(-30, 30)};
    const Vec3 start{3, rng.range(2, 8), rng.range(-3, 3)};
    const Quat q = Quat::from_axis_angle({rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)}, rng.range(0, 3));
    const Vec3 half{size, size * rng.range(0.3f, 1.6f), size * rng.range(0.3f, 1.6f)};
    World w;
    w.gravity = {};
    w.add(Body::fixed_box({thick * 0.5f, 60, 60}, {8, 5, 0}));
    Body b = trial % 2 ? Body::solid_sphere(size, 1, start) : Body::solid_box(half, 1, start, q);
    b.vel = vel;
    b.w = spin;
    w.add(b);
    return w;
}

}  // namespace

// Shots that once went through the wall with the sweep on. 418 is a plate thinner than the wall: capping
// the overlap depth did not stop it, because that depth stops growing as the plate slices in. 492 deadlocked:
// spinning at 500 rad/s it was rewound to the same pose every step and never lost its speed into the wall.
// 2014 and 666 came out of one step with hundreds of times their energy (the gyroscopic step overshooting).
TEST(shots_that_once_got_through_or_blew_up_now_bounce_off) {
    for (int trial : {418, 492, 834, 1020, 2466, 2014, 2980, 666, 608, 2402}) {
        World w = stress_shot(trial);
        const double start_energy = double(w.bodies[1].kinetic_energy());
        double peak = 0;
        for (int s = 0; s < 240; ++s) {
            w.step(kDt);
            peak = std::max(peak, double(w.bodies[1].kinetic_energy()));
        }
        const Body& b = w.bodies[1];
        CHECK(b.pos.x < 8.0f);                    // it did not get through
        CHECK(b.vel.x < 0);                       // it is on its way back, not pinned at the wall
        CHECK(peak <= start_energy * 1.001);      // and the collision never added energy
    }
}

TEST(no_shot_gains_energy_or_gets_through) {
    int through = 0, gained = 0, away = 0;
    for (int trial = 0; trial < 400; ++trial) {
        World w = stress_shot(trial);
        const double start_energy = double(w.bodies[1].kinetic_energy());
        double peak = 0;
        for (int s = 0; s < 180; ++s) {
            w.step(kDt);
            peak = std::max(peak, double(w.bodies[1].kinetic_energy()));
        }
        through += w.bodies[1].pos.x > 8.0f;
        gained += peak > start_energy * 1.001;
        away += w.bodies[1].vel.x < 0;
    }
    CHECK(through == 0);
    CHECK(gained == 0);
    CHECK(away == 400);
}
