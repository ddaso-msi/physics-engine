#include "test.hpp"
#include <phys/world.hpp>

#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <vector>

using namespace phys;

namespace {

constexpr double kG = 9.81;

Body ground(Real half_width = 20.0f) {
    Body g(Shape::make_polygon(Polygon::box(half_width, 0.5f)), {0, -0.5f}, 0, BodyType::Static);
    g.restitution = 0;
    return g;
}

Body ball(Real radius, Vec2 pos, Real restitution = 0) {
    Body b(Shape::make_circle(radius), pos, 0);
    b.restitution = restitution;
    return b;
}

Body crate(Real half, Vec2 pos, Real angle = 0) {
    Body b(Shape::make_polygon(Polygon::box(half, half)), pos, angle);
    b.restitution = 0;
    return b;
}

void run(World& w, int steps, Real dt) {
    for (int i = 0; i < steps; ++i) w.step(dt);
}

// Two circles that just overlap (well under the solver's slop), approaching along x.
World head_on(Real ra, Real rb, Real va, Real vb, Real e) {
    World w;
    w.gravity = {};
    Body a = ball(ra, {0, 0}, e), b = ball(rb, {ra + rb - 0.001f, 0}, e);
    a.vel = {va, 0};
    b.vel = {vb, 0};
    w.add(a);
    w.add(b);
    return w;
}

}  // namespace

// ---- impulses between two free bodies ------------------------------------------------------------

TEST(elastic_equal_masses_exchange_velocities) {
    World w = head_on(0.5f, 0.5f, 3, -1, 1.0f);
    w.step(1.0f / 60.0f);
    CHECK_NEAR(w.bodies[0].vel.x, -1, 1e-3);
    CHECK_NEAR(w.bodies[1].vel.x, 3, 1e-3);
}

TEST(collision_conserves_momentum_and_applies_restitution) {
    const double e = 0.3;
    World w = head_on(1.0f, 0.5f, 2, -3, static_cast<Real>(e));
    const Body &a = w.bodies[0], &b = w.bodies[1];
    const double before = a.mass * double(a.vel.x) + b.mass * double(b.vel.x);
    const double closing = double(a.vel.x) - double(b.vel.x);  // 5
    w.step(1.0f / 60.0f);
    CHECK_NEAR(a.mass * double(a.vel.x) + b.mass * double(b.vel.x), before, 1e-3);
    CHECK_NEAR(double(b.vel.x) - double(a.vel.x), e * closing, 1e-3);  // separate at e x the closing speed
}

TEST(inelastic_collision_sticks) {
    World w = head_on(0.5f, 0.5f, 4, 0, 0.0f);
    w.step(1.0f / 60.0f);
    CHECK_NEAR(w.bodies[0].vel.x, 2, 1e-3);  // equal masses share the momentum
    CHECK_NEAR(w.bodies[1].vel.x, 2, 1e-3);
}

TEST(slow_impacts_do_not_bounce) {
    // Closing at 0.5 m/s is under the restitution threshold, so even e = 1 behaves as e = 0.
    World w = head_on(0.5f, 0.5f, 0.5f, 0, 1.0f);
    w.step(1.0f / 60.0f);
    CHECK_NEAR(w.bodies[1].vel.x - w.bodies[0].vel.x, 0, 1e-3);
}

TEST(separating_contacts_are_not_pulled_back) {
    World w;
    w.add(ground());
    Body b = ball(0.5f, {0, 0.499f});  // overlapping by 1 mm
    b.vel = {0, 2};                    // but already moving away
    w.add(b);
    w.step(1.0f / 60.0f);
    CHECK_NEAR(w.bodies[1].vel.y, 2 - kG / 60.0, 1e-4);  // only gravity acted
}

TEST(static_bodies_never_move_and_static_pairs_never_collide) {
    World w;
    w.add(ground());
    w.add(Body(Shape::make_polygon(Polygon::box(1, 1)), {0, 0}, 0.3f, BodyType::Static));  // overlaps the ground
    w.add(ball(0.5f, {0, 3}));
    run(w, 240, 1.0f / 120.0f);
    CHECK(w.bodies[0].pos.x == 0 && w.bodies[0].pos.y == -0.5f);
    for (const ContactPair& c : w.contacts()) CHECK(c.a == 2 || c.b == 2);  // nothing static-vs-static
}

// ---- bouncing and resting ------------------------------------------------------------------------

TEST(bounce_follows_restitution) {
    World w;
    w.add(ground());
    const int b = w.add(ball(0.5f, {0, 3.5f}, 0.6f));  // bottom starts 3 m above the ground
    const Real dt = 1.0f / 240.0f;

    double pre = 0, post = 0, apex = 0;
    bool bounced = false;
    for (int i = 0; i < 6000; ++i) {
        const Real vy_prev = w.bodies[static_cast<size_t>(b)].vel.y;
        w.step(dt);
        const Body& ball_body = w.bodies[static_cast<size_t>(b)];
        if (!bounced && vy_prev < -1 && ball_body.vel.y > 0) {
            bounced = true;
            pre = double(vy_prev);
            post = double(ball_body.vel.y);
        }
        if (bounced) apex = std::max(apex, double(ball_body.pos.y) - 0.5);
        if (bounced && vy_prev > 0 && ball_body.vel.y <= 0) break;
    }
    CHECK(bounced);
    CHECK_NEAR(pre, -std::sqrt(2 * kG * 3), 0.1);   // free fall from 3 m
    CHECK_NEAR(post / -pre, 0.6, 0.02);             // leaves at e x arrival speed
    CHECK_NEAR(apex, 0.36 * 3, 0.06);               // and so rises to e^2 x the drop height
}

TEST(ball_comes_to_rest_on_the_ground) {
    World w;
    w.add(ground());
    const int b = w.add(ball(0.5f, {0, 2.0f}, 0.0f));
    run(w, 360, 1.0f / 120.0f);
    const Body& body = w.bodies[static_cast<size_t>(b)];
    CHECK_NEAR(body.vel.y, 0, 0.05);
    const double bottom = double(body.pos.y) - 0.5;
    CHECK(bottom > -0.015 && bottom < 0.002);  // sunk no more than about the slop
}

TEST(box_lands_flat_without_picking_up_spin) {
    World w;
    w.add(ground());
    const int b = w.add(crate(0.5f, {0, 1.5f}));
    run(w, 480, 1.0f / 120.0f);
    const Body& body = w.bodies[static_cast<size_t>(b)];
    CHECK_NEAR(body.w, 0, 0.01);
    CHECK_NEAR(body.angle, 0, 0.01);
    CHECK_NEAR(body.pos.x, 0, 0.01);
    CHECK_NEAR(body.pos.y, 0.5, 0.02);
    CHECK_NEAR(body.vel.length(), 0, 0.05);
}

// ---- friction ------------------------------------------------------------------------------------

namespace {

constexpr Real kTheta = 0.35f;  // about 20 degrees; tan = 0.365

// A big static ramp tilted by kTheta, and a box resting just above its surface, tilted to match.
World slope(Real mu) {
    World w;
    Body ramp(Shape::make_polygon(Polygon::box(20, 0.5f)), {0, 0}, kTheta, BodyType::Static);
    ramp.friction = mu;
    ramp.restitution = 0;
    w.add(ramp);
    const Vec2 normal{-std::sin(kTheta), std::cos(kTheta)};
    Body slider = crate(0.5f, normal * 1.02f, kTheta);
    slider.friction = mu;
    w.add(slider);
    return w;
}

}  // namespace

// Sliding block, Coulomb friction: a = g (sin(theta) - mu cos(theta)) down the slope.
TEST(sliding_block_accelerates_per_coulomb_friction) {
    const Real mu = 0.1f, dt = 1.0f / 120.0f;
    World w = slope(mu);
    const Vec2 downhill{-std::cos(kTheta), -std::sin(kTheta)};

    run(w, 120, dt);  // let the initial drop onto the ramp settle
    const double v1 = double(dot(w.bodies[1].vel, downhill));
    run(w, 120, dt);
    const double v2 = double(dot(w.bodies[1].vel, downhill));

    const double expected = kG * (std::sin(double(kTheta)) - double(mu) * std::cos(double(kTheta)));
    CHECK_NEAR(v2 - v1, expected, 0.05 * expected);  // one second of constant acceleration
    CHECK_NEAR(w.bodies[1].angle, kTheta, 0.02);     // slides without tipping
}

TEST(block_holds_when_friction_exceeds_tan_theta) {
    World w = slope(0.5f);  // mu = 0.5 > tan(20 deg) = 0.365
    run(w, 480, 1.0f / 120.0f);
    CHECK_NEAR(w.bodies[1].vel.length(), 0, 0.05);
}

// A circle sliding on the ground with no spin is brought to rolling by friction. For a solid disc
// (I = m r^2 / 2) the final speed is v0 / (1 + I/(m r^2)) = 2/3 v0, with w = -v/r (clockwise).
TEST(sliding_ball_settles_into_rolling) {
    World w;
    Body floor = ground();
    floor.friction = 0.5f;
    w.add(floor);
    Body b = ball(0.5f, {0, 0.5f});
    b.friction = 0.5f;
    b.vel = {6, 0};
    const int i = w.add(b);
    run(w, 180, 1.0f / 120.0f);  // slipping lasts v0/(3 mu g) = 0.41 s; this is 1.5 s
    const Body& body = w.bodies[static_cast<size_t>(i)];
    CHECK_NEAR(body.vel.x, 4, 0.1);
    CHECK_NEAR(body.w, -8, 0.3);
}

// ---- stacking and settling: where the solver's structure matters ---------------------------------

namespace {

// Four unit crates stacked with 1 mm gaps, simulated for `seconds`.
World tower(int iterations, Real seconds) {
    World w;
    w.solver.iterations = iterations;
    w.add(ground());
    for (int i = 0; i < 4; ++i) w.add(crate(0.5f, {0, 0.5f + static_cast<Real>(i) * 1.001f}));
    run(w, static_cast<int>(seconds * 120), 1.0f / 120.0f);
    return w;
}

}  // namespace

TEST(tower_of_crates_stands) {
    World w = tower(10, 5.0f);
    for (int i = 1; i <= 4; ++i) {
        const Body& b = w.bodies[static_cast<size_t>(i)];
        CHECK_NEAR(b.pos.y, 0.5 + (i - 1), 0.03);  // each crate within a few slops of its ideal height
        CHECK_NEAR(b.pos.x, 0, 0.008);             // no sideways creep
        CHECK_NEAR(b.w, 0, 0.01);
    }
}

// Each sweep fixes one contact and slightly breaks its neighbours, so a stack needs many sweeps in a
// single frame to pass weight down it... unless the solver starts each frame from the previous
// frame's impulses (warm starting). Then a single sweep per frame is enough, because the answer
// only has to be corrected a little each time rather than found from scratch.
TEST(warm_starting_lets_a_tower_stand_on_one_sweep) {
    World warm = tower(1, 5.0f);
    CHECK(warm.bodies[4].pos.y > 3.4);  // the top crate belongs at 3.5

    World cold;
    cold.solver.warm_starting = false;
    cold.solver.iterations = 1;
    cold.add(ground());
    for (int i = 0; i < 4; ++i) cold.add(crate(0.5f, {0, 0.5f + static_cast<Real>(i) * 1.001f}));
    run(cold, 600, 1.0f / 120.0f);
    CHECK(cold.bodies[4].pos.y < 3.0);  // without it the weight never reaches the bottom: the tower falls
}

// Landing on a corner and tipping flat, the box must finish truly at rest. Clamping the ACCUMULATED
// impulse (not each sweep's increment) is what lets a later sweep take back an earlier overshoot.
TEST(tilted_box_settles_flat_and_still) {
    World w;
    w.add(ground());
    w.add(crate(0.5f, {0, 2.0f}, 0.4f));
    run(w, 600, 1.0f / 120.0f);
    const Body& b = w.bodies[1];
    CHECK_NEAR(b.vel.length(), 0, 0.001);
    CHECK_NEAR(b.w, 0, 0.001);
    CHECK_NEAR(b.angle, 0, 0.01);
    CHECK_NEAR(b.pos.y, 0.5, 0.01);
}

// ---- warm starting and contact persistence -------------------------------------------------------

namespace {

ContactPair pair_with_points(int a, int b, std::initializer_list<std::uint32_t> ids, Real normal, Real tangent) {
    ContactPair p;
    p.a = a;
    p.b = b;
    for (std::uint32_t id : ids) {
        p.manifold.points[p.manifold.count].id = id;
        p.manifold.points[p.manifold.count].normal_impulse = normal;
        p.manifold.points[p.manifold.count].tangent_impulse = tangent;
        ++p.manifold.count;
    }
    return p;
}

}  // namespace

TEST(impulses_follow_the_same_body_pair_and_feature_id) {
    std::vector<ContactPair> previous = {pair_with_points(1, 2, {5}, 3.0f, 0.5f), pair_with_points(1, 3, {5}, 9.0f, 0)};
    previous[0].manifold.points[1] = previous[0].manifold.points[0];
    previous[0].manifold.points[1].id = 7;
    previous[0].manifold.points[1].normal_impulse = 4.0f;
    previous[0].manifold.points[1].tangent_impulse = -0.25f;
    previous[0].manifold.count = 2;

    // This step: pair (1,2) kept feature 7 and gained a new feature 9; pair (1,3) now touches via a
    // different feature; pair (2,3) is brand new. Stale values pre-filled to prove they get reset.
    std::vector<ContactPair> current = {pair_with_points(1, 2, {7, 9}, 99.0f, 99.0f),
                                        pair_with_points(1, 3, {6}, 99.0f, 99.0f), pair_with_points(2, 3, {5}, 0, 0)};
    transfer_impulses(previous, current);

    CHECK_NEAR(current[0].manifold.points[0].normal_impulse, 4.0, 1e-6);    // id 7 inherited
    CHECK_NEAR(current[0].manifold.points[0].tangent_impulse, -0.25, 1e-6);
    CHECK_NEAR(current[0].manifold.points[1].normal_impulse, 0, 1e-6);      // id 9 is new
    CHECK_NEAR(current[1].manifold.points[0].normal_impulse, 0, 1e-6);      // different feature: no carry-over
    CHECK_NEAR(current[2].manifold.points[0].normal_impulse, 0, 1e-6);      // different pair: no carry-over
}

TEST(resting_box_carries_the_impulse_that_cancels_gravity) {
    World w;
    w.add(ground());
    w.add(crate(0.5f, {0, 0.5f}));  // unit crate: mass 1
    const Real dt = 1.0f / 120.0f;
    run(w, 240, dt);

    CHECK(w.contacts().size() == 1);
    const Manifold& m = w.contacts()[0].manifold;
    CHECK(m.count == 2);
    // Every step gravity adds m g dt of downward velocity; the contact impulses must remove exactly that.
    double total = 0;
    for (int k = 0; k < m.count; ++k) {
        CHECK(m.points[k].normal_impulse > 0);
        total += double(m.points[k].normal_impulse);
    }
    CHECK_NEAR(total, kG * double(dt), 0.03 * kG * double(dt));
    // And they share the load equally between the two corners.
    CHECK_NEAR(m.points[0].normal_impulse, m.points[1].normal_impulse, 0.1 * double(m.points[0].normal_impulse));
}

TEST(truncate_forgets_contacts_so_new_bodies_start_fresh) {
    World w;
    w.add(ground());
    w.add(crate(0.5f, {0, 0.5f}));
    run(w, 120, 1.0f / 120.0f);
    CHECK(!w.contacts().empty());
    w.truncate(1);
    CHECK(w.bodies.size() == 1);
    CHECK(w.contacts().empty());
    w.add(crate(0.5f, {0, 1.5f}));  // reuses the old index; must not inherit the old crate's impulse
    w.step(1.0f / 120.0f);
    CHECK(w.contacts().empty());
}

// A tower of ten, only four sweeps per frame: the cold solver cannot carry the weight, the warm one can.
TEST(warm_starting_holds_a_tall_tower) {
    auto top_after = [](bool warm) {
        World w;
        w.solver.iterations = 4;
        w.solver.warm_starting = warm;
        w.add(ground());
        for (int i = 0; i < 10; ++i) w.add(crate(0.5f, {0, 0.5f + static_cast<Real>(i) * 1.001f}));
        run(w, 8 * 120, 1.0f / 120.0f);
        return double(w.bodies[10].pos.y);
    };
    CHECK(top_after(true) > 9.3);   // belongs at 9.5
    CHECK(top_after(false) < 8.0);  // collapsed
}

// Without warm starting a resting tower creeps sideways forever (the sweeps visit contacts in a fixed
// order, which biases every frame the same way). With it the tower comes to a true standstill.
TEST(resting_tower_stops_creeping) {
    World w;
    w.add(ground());
    for (int i = 0; i < 4; ++i) w.add(crate(0.5f, {8.0f, 0.5f + static_cast<Real>(i) * 1.001f}));
    run(w, 5 * 120, 1.0f / 120.0f);
    std::vector<Vec2> at5;
    for (const Body& b : w.bodies) at5.push_back(b.pos);
    run(w, 5 * 120, 1.0f / 120.0f);
    for (size_t i = 1; i < w.bodies.size(); ++i) {
        CHECK_NEAR(w.bodies[i].pos.x, at5[i].x, 0.0005);  // moved less than half a millimetre in 5 s
        CHECK_NEAR(w.bodies[i].vel.length(), 0, 0.001);
    }
}

TEST(pyramid_of_21_crates_stands) {
    World w;
    w.add(ground(30));
    const int base = 6;
    for (int row = 0; row < base; ++row)
        for (int i = 0; i < base - row; ++i)
            w.add(crate(0.5f, {(static_cast<Real>(i) - static_cast<Real>(base - row - 1) * 0.5f) * 1.02f,
                               0.5f + static_cast<Real>(row) * 1.002f}));
    std::vector<Vec2> start;
    for (const Body& b : w.bodies) start.push_back(b.pos);
    run(w, 8 * 120, 1.0f / 120.0f);
    double max_move = 0, max_speed = 0;
    for (size_t i = 1; i < w.bodies.size(); ++i) {
        max_move = std::max(max_move, double((w.bodies[i].pos - start[i]).length()));
        max_speed = std::max(max_speed, double(w.bodies[i].vel.length()));
    }
    CHECK(w.bodies.size() == 22);
    CHECK(max_move < 0.05);
    CHECK(max_speed < 0.01);
}
