#include "test.hpp"
#include <scenes.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace phys;
using namespace phys::scenes;

namespace {

constexpr double kG = 9.81;
constexpr Real kDt = 1.0f / 120.0f;

Body crate(Vec2 p, Real half = 0.5f) {
    Body b(Shape::make_polygon(Polygon::box(half, half)), p, 0);
    b.restitution = 0;
    return b;
}
Body disc(Real r, Vec2 p) {
    Body b(Shape::make_circle(r), p, 0);
    b.restitution = 0;
    return b;
}
World with_ground(Real half_width = 30) {
    World w;
    add_static_box(w, {half_width, 0.5f}, {0, -0.5f});
    return w;
}
void run(World& w, int steps) {
    for (int i = 0; i < steps; ++i) w.step(kDt);
}
// Steps until `pred` or `max_steps`; returns the step count at which it became true, or -1.
template <class Pred>
int run_until(World& w, int max_steps, Pred pred) {
    for (int i = 0; i < max_steps; ++i) {
        w.step(kDt);
        if (pred()) return i + 1;
    }
    return -1;
}

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
};

// Bodies 1..n of a pyramid of crates (body 0 is the ground).
void add_pyramid(World& w, int base) {
    for (int row = 0; row < base; ++row)
        for (int i = 0; i < base - row; ++i)
            w.add(crate({(static_cast<Real>(i) - static_cast<Real>(base - row - 1) * 0.5f) * 1.02f,
                         0.5f + static_cast<Real>(row) * 1.002f}));
}

}  // namespace

// ---- sleeping ------------------------------------------------------------------------------------

TEST(a_resting_box_falls_asleep_and_then_does_not_move_at_all) {
    World w = with_ground();
    const int c = w.add(crate({0, 1.5f}));
    run(w, 360);  // dropped, bounced about, settled, still for half a second
    CHECK(!w.bodies[static_cast<size_t>(c)].awake);
    CHECK(w.stats().awake_bodies == 0);
    CHECK(w.bodies[static_cast<size_t>(c)].vel.length_sq() == 0);

    const Vec2 p = w.bodies[static_cast<size_t>(c)].pos;
    const Real angle = w.bodies[static_cast<size_t>(c)].angle;
    run(w, 240);
    CHECK(w.bodies[static_cast<size_t>(c)].pos.x == p.x && w.bodies[static_cast<size_t>(c)].pos.y == p.y);  // bit for bit
    CHECK(w.bodies[static_cast<size_t>(c)].angle == angle);

    // Its contact with the ground is remembered (dormant) along with the impulse that held it up.
    CHECK(!w.contacts().empty());
    double total = 0;
    for (const ContactPair& cp : w.contacts()) {
        CHECK(cp.dormant);
        for (int k = 0; k < cp.manifold.count; ++k) total += double(cp.manifold.points[k].normal_impulse);
    }
    CHECK_NEAR(total, kG * double(kDt), 0.05 * kG * double(kDt));  // m g dt, unit mass
    CHECK(w.stats().contacts == 0);  // nothing was solved
}

TEST(bodies_that_keep_moving_do_not_sleep) {
    World w;  // a pendulum never comes to rest in four seconds
    w.add(disc(0.2f, {-1.0f, 5}));
    w.add_joint(Joint::revolute(w.bodies, -1, 0, {0, 5}));
    run(w, 480);
    CHECK(w.bodies[0].awake);
}

TEST(sleep_thresholds_separate_slow_from_still) {
    auto asleep_after_1s = [](Real speed) {
        World w;
        w.gravity = {};
        Body b = disc(0.2f, {0, 0});
        b.vel = {speed, 0};
        w.add(b);
        run(w, 120);
        return !w.bodies[0].awake;
    };
    CHECK(!asleep_after_1s(0.02f));  // above the 0.01 m/s threshold: drifting, not resting
    CHECK(asleep_after_1s(0.005f));  // below it: asleep after half a second
}

TEST(an_island_sleeps_and_wakes_as_one) {
    World w = with_ground();
    for (int i = 0; i < 4; ++i) w.add(crate({0, 0.5f + static_cast<Real>(i) * 1.001f}));
    // The tower must go to sleep all at once, never crate by crate.
    int steps = 0;
    for (; steps < 20 * 120; ++steps) {
        w.step(kDt);
        int asleep = 0;
        for (int i = 1; i <= 4; ++i) asleep += !w.bodies[static_cast<size_t>(i)].awake;
        CHECK(asleep == 0 || asleep == 4);
        if (asleep == 4) break;
    }
    CHECK(steps < 20 * 120);

    // A ball dropped on the top crate wakes the whole tower in the same step, before the ball even
    // reaches the lower crates.
    const int ball = w.add(disc(0.3f, {0, 7}));
    const int woke = run_until(w, 240, [&] { return w.bodies[1].awake || w.bodies[2].awake || w.bodies[3].awake || w.bodies[4].awake; });
    CHECK(woke > 0);
    for (int i = 1; i <= 4; ++i) CHECK(w.bodies[static_cast<size_t>(i)].awake);
    CHECK(w.bodies[static_cast<size_t>(ball)].awake);
}

TEST(a_sleeping_pyramid_is_not_disturbed_by_a_poke) {
    World w = with_ground();
    add_pyramid(w, 5);
    run(w, 6 * 120);
    CHECK(w.stats().awake_bodies == 0);

    std::vector<Vec2> before;
    for (const Body& b : w.bodies) before.push_back(b.pos);
    const size_t top = w.bodies.size() - 1;
    w.bodies[top].apply_impulse_at({0.05f, 0}, w.bodies[top].pos);  // a tap on the top crate
    w.step(kDt);
    CHECK(w.bodies[top].awake);  // the tap woke it (and with it the pile)
    run(w, 120);
    double worst_below = 0;
    for (size_t i = 1; i < top; ++i) worst_below = std::max(worst_below, double((w.bodies[i].pos - before[i]).length()));
    CHECK(worst_below < 0.02);  // the pile woke up warm started: no jump, no sag
}

TEST(setting_velocity_or_force_on_a_sleeper_wakes_it) {
    World w = with_ground();
    const int a = w.add(crate({-3, 0.5f}));
    const int b = w.add(crate({3, 0.5f}));
    run(w, 240);
    CHECK(!w.bodies[static_cast<size_t>(a)].awake && !w.bodies[static_cast<size_t>(b)].awake);

    w.bodies[static_cast<size_t>(a)].vel = {2, 0};
    w.bodies[static_cast<size_t>(b)].apply_force({200, 0});
    run(w, 30);
    CHECK(w.bodies[static_cast<size_t>(a)].pos.x > -3 + 0.1);
    CHECK(w.bodies[static_cast<size_t>(b)].pos.x > 3 + 0.1);
}

TEST(a_body_moved_by_hand_needs_wake_to_respond) {
    World w;
    const int b = w.add(disc(0.3f, {0, 2}));
    w.gravity = {};
    run(w, 120);
    CHECK(!w.bodies[static_cast<size_t>(b)].awake);
    w.gravity = {0, static_cast<Real>(-kG)};
    w.bodies[static_cast<size_t>(b)].pos = {0, 10};  // teleported, nobody told the world
    run(w, 60);
    CHECK(w.bodies[static_cast<size_t>(b)].pos.y == 10);  // still frozen
    w.wake(b);
    run(w, 60);
    CHECK(w.bodies[static_cast<size_t>(b)].pos.y < 9.9f);  // now it falls
}

TEST(a_mouse_joint_wakes_a_sleeper_and_drags_it) {
    World w = with_ground();
    const int c = w.add(crate({0, 0.5f}));
    run(w, 240);
    CHECK(!w.bodies[static_cast<size_t>(c)].awake);
    Joint grab = Joint::mouse(w.bodies, c, w.bodies[static_cast<size_t>(c)].pos);
    grab.target = {2, 0.5f};
    w.add_joint(grab);
    run(w, 180);
    CHECK(w.bodies[static_cast<size_t>(c)].pos.x > 1.0f);
}

TEST(a_motor_told_to_turn_wakes_a_sleeping_mechanism) {
    World w;
    w.gravity = {};
    w.add(disc(0.5f, {0, 3}));
    Joint axle = Joint::revolute(w.bodies, -1, 0, {0, 3});
    axle.enable_motor = true;
    axle.max_motor = 100;
    axle.motor_speed = 0;
    w.add_joint(axle);
    run(w, 120);
    CHECK(!w.bodies[0].awake);
    w.joints[0].motor_speed = 2;
    run(w, 60);
    CHECK_NEAR(w.bodies[0].w, 2, 0.05);
}

TEST(jointed_bodies_share_an_island) {
    World w = with_ground();
    const int a = w.add(disc(0.3f, {-0.75f, 0.3f}));
    const int b = w.add(disc(0.3f, {0.75f, 0.3f}));
    w.add_joint(Joint::distance(w.bodies, a, b, w.bodies[static_cast<size_t>(a)].pos, w.bodies[static_cast<size_t>(b)].pos));
    run(w, 240);
    CHECK(!w.bodies[static_cast<size_t>(a)].awake && !w.bodies[static_cast<size_t>(b)].awake);
    w.bodies[static_cast<size_t>(a)].apply_impulse_at({0.2f, 0}, w.bodies[static_cast<size_t>(a)].pos);
    w.step(kDt);
    CHECK(w.bodies[static_cast<size_t>(b)].awake);  // pulled along by the rod
    run(w, 20);
    CHECK(w.bodies[static_cast<size_t>(b)].pos.x > 0.75f + 0.005f);
}

TEST(islands_are_counted_and_static_bodies_do_not_join_them) {
    World apart = with_ground();
    apart.add(crate({-3, 0.5f}));
    apart.add(crate({3, 0.5f}));
    apart.step(kDt);
    CHECK(apart.stats().islands == 2);  // both sit on the same floor but are not linked through it

    World touching = with_ground();
    touching.add(crate({0, 0.5f}));
    touching.add(crate({0.99f, 0.5f}));  // overlapping a little
    touching.step(kDt);
    CHECK(touching.stats().islands == 1);
}

TEST(a_body_that_may_not_sleep_keeps_its_whole_island_awake) {
    World w = with_ground();
    w.add(crate({0, 0.5f}));
    Body top = crate({0, 1.501f});
    top.allow_sleep = false;
    w.add(top);
    run(w, 600);
    CHECK(w.bodies[1].awake && w.bodies[2].awake);
    w.allow_sleep = false;
    run(w, 10);
    CHECK(w.stats().awake_bodies == 2);
}

TEST(sleeping_saves_work_without_changing_the_result) {
    World sleepy = with_ground(), awake = with_ground();
    awake.allow_sleep = false;
    add_pyramid(sleepy, 6);
    add_pyramid(awake, 6);
    run(sleepy, 8 * 120);
    run(awake, 8 * 120);
    CHECK(sleepy.stats().awake_bodies == 0);
    CHECK(awake.stats().awake_bodies == 21);
    CHECK(sleepy.stats().contacts == 0);   // nothing left to solve...
    CHECK(awake.stats().contacts > 20);    // ...where the awake world still solves them all
    double worst = 0;
    for (size_t i = 0; i < sleepy.bodies.size(); ++i)
        worst = std::max(worst, double((sleepy.bodies[i].pos - awake.bodies[i].pos).length()));
    CHECK(worst < 0.02);
}

TEST(truncating_a_world_with_sleepers_is_safe) {
    World w = with_ground();
    add_pyramid(w, 4);
    run(w, 6 * 120);
    CHECK(w.stats().awake_bodies == 0);
    w.truncate(1);
    CHECK(w.contacts().empty());
    w.add(crate({0, 3}));
    run(w, 240);
    CHECK_NEAR(w.bodies[1].pos.y, 0.5, 0.02);
}

// ---- continuous collision ------------------------------------------------------------------------

namespace {

// A tall, thin static wall at x = 8: 4 cm thick.
void add_thin_wall(World& w) { add_static_box(w, {0.02f, 20}, {8, 5}); }

}  // namespace

TEST(a_fast_ball_stops_at_a_thin_wall_only_with_continuous_collision) {
    auto final_x = [](bool continuous) {
        World w;
        w.gravity = {};
        w.continuous = continuous;
        add_thin_wall(w);
        Body bullet = disc(0.15f, {2, 5});
        bullet.restitution = 0.5f;
        bullet.vel = {60, 0};  // 0.5 m per step: more than the wall and ball together
        w.add(bullet);
        run(w, 120);
        return std::pair<Real, Real>(w.bodies[1].pos.x, w.bodies[1].vel.x);
    };
    const auto without = final_x(false);
    CHECK(without.first > 8.5f);     // sailed straight through
    const auto with = final_x(true);
    CHECK(with.first < 8.0f);        // stopped on this side...
    CHECK(with.second < 0);          // ...and bounced back
    CHECK_NEAR(with.second, -30, 3); // at e x the arrival speed
}

TEST(continuous_collision_stops_a_spinning_box) {
    World w;
    w.gravity = {};
    add_thin_wall(w);
    Body bullet(Shape::make_polygon(Polygon::box(0.3f, 0.05f)), {2, 5}, 0.4f);
    bullet.vel = {80, 3};
    bullet.w = 25;
    w.add(bullet);
    run(w, 60);
    CHECK(w.bodies[1].pos.x < 8.0f);
}

TEST(continuous_collision_works_for_any_shot) {
    int through_with = 0, through_without = 0;
    for (int trial = 0; trial < 150; ++trial) {
        Lcg rng;
        rng.s = static_cast<std::uint32_t>(trial) * 7919u + 3;
        const Real speed = rng.range(40, 120), vy = rng.range(-25, 25), spin = rng.range(-12, 12), y = rng.range(2, 8);
        const bool box = trial % 2 == 0;
        for (int mode = 0; mode < 2; ++mode) {
            World w;
            w.gravity = {};
            w.continuous = mode == 1;
            add_thin_wall(w);
            Body bullet = box ? Body(Shape::make_polygon(Polygon::box(0.12f, 0.2f)), {3, y}, rng.range(-3, 3))
                              : Body(Shape::make_circle(0.12f), {3, y}, 0);
            bullet.vel = {speed, vy};
            bullet.w = spin;
            w.add(bullet);
            run(w, 60);
            const bool through = w.bodies[1].pos.x > 8.0f;
            (mode == 1 ? through_with : through_without) += through;
        }
    }
    CHECK(through_with == 0);       // nothing ever gets through
    CHECK(through_without > 90);    // and without it most shots do: the test is not vacuous
}

namespace {

// One of a family of much nastier shots: plates as thin as 1:10, walls from 1 to 20 cm, up to 200 m/s and
// 80 rad/s. `kind` 0 is a disc, 1 a box, 2 a triangle whose centre of mass is off to one side.
World harsh_shot(int trial) {
    Lcg rng;
    rng.s = static_cast<std::uint32_t>(trial) * 104729u + 11;
    const Real speed = rng.range(20, 200), vy = rng.range(-40, 40), spin = rng.range(-80, 80), y = rng.range(2, 8);
    const Real thick = rng.range(0.01f, 0.2f), size = rng.range(0.06f, 0.4f), aspect = rng.range(0.1f, 1.6f);
    const Real angle = rng.range(-3, 3);
    World w;
    w.gravity = {};
    w.allow_sleep = false;
    add_static_box(w, {thick * 0.5f, 60}, {8, 5});
    const int kind = trial % 3;
    const Vec2 corners[3] = {{0, 0}, {size * 2, 0}, {size * 0.5f, size * 2 * aspect}};
    Body bullet = kind == 0   ? Body(Shape::make_circle(size), {3, y}, 0)
                  : kind == 1 ? Body(Shape::make_polygon(Polygon::box(size, size * aspect)), {3, y}, angle)
                              : Body(Shape::make_polygon(*Polygon::from_points(corners)), {3, y}, angle);
    bullet.vel = {speed, vy};
    bullet.w = spin;
    w.add(bullet);
    return w;
}

}  // namespace

// Three ways to get this wrong: let the body through, hold it against the wall for ever with its speed
// intact, or hand it energy. An earlier version of the sweep managed the second, one shot in ten.
TEST(continuous_collision_neither_leaks_nor_traps_nor_adds_energy) {
    int through = 0, trapped = 0, gained = 0;
    for (int trial = 0; trial < 600; ++trial) {
        World w = harsh_shot(trial);
        const Real e0 = w.bodies[1].kinetic_energy();
        Real peak = 0;
        for (int s = 0; s < 240; ++s) {
            w.step(1.0f / 120.0f);
            peak = std::max(peak, w.bodies[1].kinetic_energy());
        }
        const Body& b = w.bodies[1];
        through += b.pos.x > 8.3f;
        trapped += b.pos.x < 8.0f && b.vel.x >= 0.5f;  // still this side, two seconds on, still heading in
        gained += peak > e0 * 1.001f;
    }
    CHECK(through == 0);
    CHECK(trapped == 0);
    CHECK(gained == 0);
}

// The shot that exposed the trap. An off-centre hit turns the box's speed into spin: 407 rad/s, more than
// half a revolution per step. The sweep used to work out the turn from the two stored angles, which are
// wrapped, so it swept the box round the wrong way and put it back in the same pose every step.
TEST(a_box_set_spinning_by_its_impact_still_comes_back_off_the_wall) {
    World w = harsh_shot(1);
    for (int i = 0; i < 240; ++i) w.step(1.0f / 120.0f);
    CHECK(w.bodies[1].pos.x < 8.0f);
    CHECK(w.bodies[1].vel.x < 0);
}

TEST(a_fast_fall_onto_a_thin_floor_does_not_tunnel) {
    World w;
    Body floor(Shape::make_polygon(Polygon::box(20, 0.02f)), {0, 0}, 0, BodyType::Static);
    floor.restitution = 0;
    w.add(floor);
    Body ball = disc(0.2f, {0, 6});
    ball.vel = {0, -150};
    w.add(ball);
    run(w, 120);
    CHECK(w.bodies[1].pos.y > 0.1f);
}

TEST(continuous_collision_leaves_slow_bodies_alone) {
    World a = with_ground(), b = with_ground();
    b.continuous = false;
    add_pyramid(a, 5);
    add_pyramid(b, 5);
    for (int i = 0; i < 3; ++i) {
        a.add(disc(0.3f, {-2 + static_cast<Real>(i) * 2, 6}));
        b.add(disc(0.3f, {-2 + static_cast<Real>(i) * 2, 6}));
    }
    run(a, 240);
    run(b, 240);
    CHECK(a.stats().ccd_hits == 0);
    for (size_t i = 0; i < a.bodies.size(); ++i)
        CHECK(a.bodies[i].pos.x == b.bodies[i].pos.x && a.bodies[i].pos.y == b.bodies[i].pos.y);  // bit for bit
}

TEST(continuous_collision_respects_collision_filters) {
    World w;
    w.gravity = {};
    add_thin_wall(w);
    Body ghost = disc(0.15f, {2, 5});
    ghost.category = 0b10;
    ghost.mask = 0;  // collides with nothing
    ghost.vel = {60, 0};
    w.add(ghost);
    run(w, 60);
    CHECK(w.bodies[1].pos.x > 8.5f);
    CHECK(w.stats().ccd_hits == 0);
}

// ---- the gallery ---------------------------------------------------------------------------------

TEST(newtons_cradle_passes_the_swing_to_the_far_ball) {
    World w;
    const Cradle c = build_cradle(w, {8, 9.6f});
    lift_cradle_ball(w, c, 0, -0.7f);
    const double rest = 9.6 - 4.0;
    const double gain = double(w.bodies[static_cast<size_t>(c.balls[0])].pos.y) - rest;
    std::vector<double> peak(5, -1e9);
    for (int s = 0; s < 3 * 240; ++s) {
        w.step(1.0f / 240.0f);
        for (size_t i = 0; i < 5; ++i) peak[i] = std::max(peak[i], double(w.bodies[static_cast<size_t>(c.balls[i])].pos.y));
    }
    CHECK(gain > 0.9);
    CHECK(peak[4] - rest > 0.95 * gain);   // the far ball swings out almost as high as the first was lifted
    for (size_t i = 1; i <= 3; ++i) CHECK(peak[i] - rest < 0.02);  // the middle balls hardly stir
}

TEST(two_lifted_balls_send_two_out_the_far_side) {
    World w;
    const Cradle c = build_cradle(w, {8, 9.6f});
    lift_cradle_ball(w, c, 0, -0.6f);
    lift_cradle_ball(w, c, 1, -0.6f);
    const double rest = 9.6 - 4.0, gain = double(w.bodies[static_cast<size_t>(c.balls[0])].pos.y) - rest;
    std::vector<double> peak(5, -1e9);
    for (int s = 0; s < 3 * 240; ++s) {
        w.step(1.0f / 240.0f);
        for (size_t i = 0; i < 5; ++i) peak[i] = std::max(peak[i], double(w.bodies[static_cast<size_t>(c.balls[i])].pos.y));
    }
    CHECK(peak[3] - rest > 0.9 * gain);
    CHECK(peak[4] - rest > 0.9 * gain);
    CHECK(peak[2] - rest < 0.05);  // the middle one stays put
}

TEST(a_ragdoll_falls_folds_within_its_limits_and_comes_to_rest) {
    World w = with_ground(20);
    const Ragdoll d = build_ragdoll(w, {0, kRagdollTorsoHeight + 3.0f}, 1);
    double worst_gap = 0, lowest = 1e9;
    double worst_limit_violation = 0;
    int asleep_at = -1;
    for (int s = 0; s < 12 * 120; ++s) {
        w.step(kDt);
        for (const Joint& j : w.joints) {
            worst_gap = std::max(worst_gap, double((j.world_anchor_b(w.bodies) - j.world_anchor_a(w.bodies)).length()));
            const double angle = double(w.bodies[static_cast<size_t>(j.b)].angle) - double(w.bodies[static_cast<size_t>(j.a)].angle) - double(j.reference_angle);
            worst_limit_violation = std::max({worst_limit_violation, double(j.lower) - angle, angle - double(j.upper)});
        }
        for (const Body& b : w.bodies)
            if (b.type == BodyType::Dynamic) lowest = std::min(lowest, double(b.pos.y));
        if (asleep_at < 0 && w.stats().awake_bodies == 0) asleep_at = s;
    }
    CHECK(w.joints.size() == 9 && d.hinges.size() == 9);
    CHECK(worst_gap < 0.08);               // limbs stay attached through the impact
    CHECK(worst_limit_violation < 0.05);   // elbows, knees, neck all stay inside their limits
    CHECK(lowest > 0.05);                  // nothing sank into the floor
    CHECK(asleep_at > 0 && asleep_at < 8 * 120);  // it lies still and goes to sleep
    // It ended up lying down: the torso (which stood 1.74 up, 1.1 tall) is low and tipped over. The head
    // does not reach the floor, because the neck's limit holds it up off a torso propped on its arms.
    const Body& torso = w.bodies[static_cast<size_t>(d.torso)];
    CHECK(torso.pos.y < 0.6f);
    CHECK(std::fabs(torso.angle) > 1.0f);
    CHECK(w.bodies[static_cast<size_t>(d.head)].pos.y < 0.9f);
}

TEST(a_car_settles_on_its_springs_drives_and_brakes) {
    World w = with_ground(80);
    const Car car = build_car(w, {0, 1.1f}, 1);
    run(w, 4 * 120);
    const Body& chassis = w.bodies[static_cast<size_t>(car.chassis)];
    CHECK(chassis.pos.y > 0.85f && chassis.pos.y < 1.02f);  // sagged on the springs, not bottomed out (1.07 unloaded, 0.82 at the stops)
    CHECK_NEAR(chassis.angle, 0, 0.01);
    CHECK(!chassis.awake);  // the parked car sleeps, so the throttle below has to wake it

    set_car_throttle(w, car, 12.0f, 25.0f);
    const double x0 = double(chassis.pos.x);
    double worst_tilt = 0, worst_axle = 0;
    for (int s = 0; s < 4 * 120; ++s) {
        w.step(kDt);
        worst_tilt = std::max(worst_tilt, double(std::fabs(chassis.angle)));
        for (const Joint& j : w.joints)
            if (j.type == JointType::Revolute) worst_axle = std::max(worst_axle, double((j.world_anchor_b(w.bodies) - j.world_anchor_a(w.bodies)).length()));
    }
    CHECK(double(chassis.pos.x) - x0 > 15.0);        // about 20 m in 4 s
    CHECK_NEAR(chassis.vel.x, 12 * 0.45, 0.4);       // wheel speed x radius
    CHECK(worst_tilt < 0.4);
    CHECK(worst_axle < 0.005);

    set_car_throttle(w, car, 0.0f, 40.0f);  // brake
    const int stopped = run_until(w, 3 * 120, [&] { return std::fabs(chassis.vel.x) < 0.05f; });
    CHECK(stopped > 0);
    run(w, 8 * 120);
    CHECK(!chassis.awake);  // and the stopped car goes back to sleep
}

// A long rod, pinned at its centre and spinning at 40 rad/s, sweeps 0.33 rad per step: its tip moves 33 cm
// a step even though the rod's centre does not move at all. A small post sits on the arc between two steps' positions, so at no
// step boundary does the rod touch it; only the sweep in between hits it. Rotation must count as motion.
TEST(a_fast_spinner_does_not_swing_through_a_post) {
    auto furthest_swing = [](bool continuous) {
        World w;
        w.gravity = {};
        w.continuous = continuous;
        const Real post_angle = 0.5f, reach = 0.9f;
        add_static_box(w, {0.04f, 0.04f}, {reach * std::cos(post_angle), reach * std::sin(post_angle)});
        Body rod(Shape::make_polygon(Polygon::box(1.0f, 0.05f)), {0, 0}, 0);
        rod.w = 40;
        rod.restitution = 0;
        w.add(rod);
        w.add_joint(Joint::revolute(w.bodies, -1, 1, {0, 0}));  // pinned at its centre: a turnstile arm that can only turn
        Real furthest = 0;
        for (int s = 0; s < 36; ++s) {
            w.step(kDt);
            furthest = std::max(furthest, w.bodies[1].angle);  // the spin stays positive until the rod is stopped
        }
        return furthest;
    };
    CHECK(furthest_swing(false) > 0.9f);  // without the sweep it goes straight through the post
    CHECK(furthest_swing(true) < 0.55f);  // with it, it is stopped on this side
}

TEST(only_fast_bodies_are_swept) {
    World pile = with_ground();
    pile.allow_sleep = false;  // keep every crate live, so "nothing swept" cannot just mean "everything asleep"
    add_pyramid(pile, 5);
    run(pile, 240);
    CHECK(pile.stats().awake_bodies == 15);
    CHECK(pile.stats().ccd_swept == 0);  // a settled pile never triggers the sweep

    World shot;
    shot.gravity = {};
    add_static_box(shot, {0.02f, 20}, {8, 5});
    Body bullet = disc(0.15f, {2, 5});
    bullet.vel = {60, 0};
    shot.add(bullet);
    shot.step(kDt);
    CHECK(shot.stats().ccd_swept == 1);
}
