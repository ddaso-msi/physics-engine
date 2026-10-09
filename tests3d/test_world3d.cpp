#include "test.hpp"
#include <phys3d/world.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

using namespace phys3d;

namespace {

constexpr double kG = 9.81;
constexpr Real kDt = 1.0f / 120.0f;

#define CHECK_VEC(a, b, eps)                         \
    do {                                             \
        const Vec3 va_ = (a), vb_ = (b);             \
        CHECK_NEAR(va_.x, vb_.x, eps);               \
        CHECK_NEAR(va_.y, vb_.y, eps);               \
        CHECK_NEAR(va_.z, vb_.z, eps);               \
    } while (0)

Body floor_box(Real friction = 0.5f) {
    Body f = Body::fixed_box({30, 0.5f, 30}, {0, -0.5f, 0});  // top face at y = 0
    f.restitution = 0;
    f.friction = friction;
    return f;
}
Body crate(Vec3 p, Real half = 0.5f) {
    Body b = Body::solid_box({half, half, half}, 1, p);
    b.restitution = 0;
    return b;
}
Body ball(Real r, Vec3 p, Real restitution = 0) {
    Body b = Body::solid_sphere(r, 1, p);
    b.restitution = restitution;
    return b;
}
void run(World& w, int steps, Real dt = kDt) {
    for (int i = 0; i < steps; ++i) w.step(dt);
}

// Two spheres barely overlapping (under the slop), closing along a line through the given direction.
World head_on(Real ra, Real rb, Real va, Real vb, Real e, Vec3 along = {1, 0, 0}) {
    World w;
    w.gravity = {};
    const Vec3 d = along.normalized();
    Body a = ball(ra, {0, 0, 0}, e), b = ball(rb, d * (ra + rb - 0.001f), e);
    a.vel = d * va;
    b.vel = d * vb;
    w.add(a);
    w.add(b);
    return w;
}

// How far the bodies after index 0 have moved from `start`, and how fast anything is still going.
struct Settle {
    double moved = 0, speed = 0;
};
Settle settle(const World& w, const std::vector<Vec3>& start) {
    Settle s;
    for (size_t i = 1; i < w.bodies.size(); ++i) {
        s.moved = std::max(s.moved, double((w.bodies[i].pos - start[i]).length()));
        s.speed = std::max(s.speed, double(w.bodies[i].vel.length()) + double(w.bodies[i].w.length()));
    }
    return s;
}
std::vector<Vec3> positions(const World& w) {
    std::vector<Vec3> p;
    for (const Body& b : w.bodies) p.push_back(b.pos);
    return p;
}

}  // namespace

// ---- impulses between two free bodies ------------------------------------------------------------

TEST(elastic_equal_spheres_exchange_velocities_along_any_line) {
    for (Vec3 dir : {Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{1, 2, -3}}) {
        World w = head_on(0.5f, 0.5f, 3, -1, 1.0f, dir);
        w.step(kDt);
        const Vec3 d = dir.normalized();
        CHECK_VEC(w.bodies[0].vel, d * -1.0f, 2e-3);
        CHECK_VEC(w.bodies[1].vel, d * 3.0f, 2e-3);
    }
}

TEST(collision_conserves_momentum_and_applies_restitution_in_3d) {
    const double e = 0.3;
    World w = head_on(1.0f, 0.5f, 2, -3, static_cast<Real>(e), {1, 1, 0});
    const Body &a = w.bodies[0], &b = w.bodies[1];
    const Vec3 before = a.vel * a.mass + b.vel * b.mass;
    w.step(kDt);
    CHECK_VEC(a.vel * a.mass + b.vel * b.mass, before, 2e-3);
    const Vec3 d = Vec3{1, 1, 0}.normalized();
    CHECK_NEAR(dot(b.vel - a.vel, d), e * 5, 2e-3);  // separate at e x the closing speed of 5
}

TEST(inelastic_spheres_stick_and_slow_ones_do_not_bounce) {
    World sticky = head_on(0.5f, 0.5f, 4, 0, 0.0f);
    sticky.step(kDt);
    CHECK_NEAR(sticky.bodies[0].vel.x, 2, 2e-3);
    CHECK_NEAR(sticky.bodies[1].vel.x, 2, 2e-3);

    World slow = head_on(0.5f, 0.5f, 0.5f, 0, 1.0f);  // closing at 0.5 m/s, under the 1 m/s threshold
    slow.step(kDt);
    CHECK_NEAR(slow.bodies[1].vel.x - slow.bodies[0].vel.x, 0, 2e-3);
}

TEST(contacts_push_but_never_pull_in_3d) {
    World w;
    w.add(floor_box());
    Body b = ball(0.5f, {0, 0.499f, 0});  // 1 mm into the floor...
    b.vel = {0, 2, 0};                    // ...but already leaving
    w.add(b);
    w.step(1.0f / 60.0f);
    CHECK_NEAR(w.bodies[1].vel.y, 2 - kG / 60.0, 1e-4);  // only gravity acted
}

TEST(static_bodies_never_move_and_static_pairs_make_no_contacts) {
    World w;
    w.add(floor_box());
    w.add(Body::fixed_box({1, 1, 1}, {0, 0.2f, 0}));  // overlapping the floor
    w.add(ball(0.5f, {5, 3, 0}));
    run(w, 240);
    CHECK_VEC(w.bodies[0].pos, Vec3(0, -0.5f, 0), 0);
    CHECK_VEC(w.bodies[1].pos, Vec3(0, 0.2f, 0), 0);
    for (const ContactPair& c : w.contacts()) CHECK(c.a == 2 || c.b == 2);
}

// ---- bouncing and resting ------------------------------------------------------------------------

TEST(a_ball_bounces_to_e_squared_of_its_drop_height) {
    World w;
    w.add(floor_box());
    const int b = w.add(ball(0.5f, {0, 3.5f, 0}, 0.6f));  // bottom 3 m above the floor
    const Real dt = 1.0f / 240.0f;
    double pre = 0, post = 0, apex = 0;
    bool bounced = false;
    for (int i = 0; i < 6000; ++i) {
        const Real vy_prev = w.bodies[static_cast<size_t>(b)].vel.y;
        w.step(dt);
        const Body& body = w.bodies[static_cast<size_t>(b)];
        if (!bounced && vy_prev < -1 && body.vel.y > 0) {
            bounced = true;
            pre = double(vy_prev);
            post = double(body.vel.y);
        }
        if (bounced) apex = std::max(apex, double(body.pos.y) - 0.5);
        if (bounced && vy_prev > 0 && body.vel.y <= 0) break;
    }
    CHECK(bounced);
    CHECK_NEAR(pre, -std::sqrt(2 * kG * 3), 0.1);
    CHECK_NEAR(post / -pre, 0.6, 0.02);
    CHECK_NEAR(apex, 0.36 * 3, 0.06);
}

TEST(a_crate_rests_on_four_points_that_carry_its_weight) {
    World w;
    w.add(floor_box());
    w.add(crate({0, 0.5f, 0}));  // unit cube, mass 1
    run(w, 240);
    CHECK(w.contacts().size() == 1);
    if (w.contacts().empty()) return;  // already failed; do not index into nothing
    const Manifold& m = w.contacts()[0].manifold;
    CHECK(m.count == 4);
    double total = 0;
    for (int k = 0; k < m.count; ++k) {
        CHECK(m.points[k].normal_impulse > 0);
        total += double(m.points[k].normal_impulse);
        // Each corner carries about a quarter (they need not be exactly equal: any split with the right
        // total and no net torque holds the crate up).
        CHECK_NEAR(m.points[k].normal_impulse, kG * double(kDt) / 4, 0.1 * kG * double(kDt) / 4);
    }
    CHECK_NEAR(total, kG * double(kDt), 0.02 * kG * double(kDt));  // m g dt every step
    const Body& c = w.bodies[1];
    CHECK(c.pos.y < 0.5f && c.pos.y > 0.49f);  // sunk by no more than about the slop
    CHECK_NEAR(double(c.vel.length()) + double(c.w.length()), 0, 1e-3);
}

TEST(a_crate_dropped_on_a_corner_ends_up_flat_and_still) {
    World w;
    w.add(floor_box());
    Body b = Body::solid_box({0.5f, 0.5f, 0.5f}, 1, {0, 2, 0}, Quat::from_axis_angle({1, 0.3f, 0.5f}, 0.6f));
    b.restitution = 0;
    w.add(b);
    run(w, 720);
    const Body& c = w.bodies[1];
    CHECK_NEAR(c.pos.y, 0.5, 0.01);
    CHECK_NEAR(double(c.vel.length()) + double(c.w.length()), 0, 1e-3);
    // One of its own axes points straight up: it is lying on a face.
    const Real up = std::max({std::fabs(rotate(c.q, {1, 0, 0}).y), std::fabs(rotate(c.q, {0, 1, 0}).y), std::fabs(rotate(c.q, {0, 0, 1}).y)});
    CHECK_NEAR(up, 1, 1e-3);
}

// ---- friction ------------------------------------------------------------------------------------

namespace {

constexpr Real kTheta = 0.35f;  // about 20 degrees; tan = 0.365

// A big ramp tilted by kTheta about z, and a crate resting just above it, tilted to match.
World slope(Real mu) {
    World w;
    Body ramp = Body::fixed_box({20, 0.5f, 20}, {0, 0, 0}, Quat::from_axis_angle({0, 0, 1}, kTheta));
    ramp.friction = mu;
    ramp.restitution = 0;
    w.add(ramp);
    const Vec3 normal{-std::sin(kTheta), std::cos(kTheta), 0};
    Body slider = Body::solid_box({0.5f, 0.5f, 0.5f}, 1, normal * 1.02f, Quat::from_axis_angle({0, 0, 1}, kTheta));
    slider.friction = mu;
    slider.restitution = 0;
    w.add(slider);
    return w;
}
const Vec3 kDownhill{-std::cos(kTheta), -std::sin(kTheta), 0};

}  // namespace

TEST(a_crate_slides_down_a_slope_per_coulomb_friction) {
    World w = slope(0.1f);
    run(w, 120);
    const double v1 = double(dot(w.bodies[1].vel, kDownhill));
    run(w, 120);
    const double v2 = double(dot(w.bodies[1].vel, kDownhill));
    const double expected = kG * (std::sin(double(kTheta)) - 0.1 * std::cos(double(kTheta)));
    CHECK_NEAR(v2 - v1, expected, 0.02 * expected);  // one second of g (sin - mu cos)
    CHECK_NEAR(w.bodies[1].vel.z, 0, 1e-3);          // straight down the slope, no sideways drift
    CHECK_NEAR(w.bodies[1].w.length(), 0, 1e-3);     // and without tumbling
}

TEST(a_crate_holds_on_a_slope_when_friction_exceeds_tan_theta) {
    World w = slope(0.5f);
    run(w, 480);
    CHECK_NEAR(w.bodies[1].vel.length(), 0, 0.02);
}

// Friction must be the same in every direction. A slab sliding along the diagonal of the two tangent axes
// decelerates at mu g like any other; clamping the two tangent rows separately would brake it 41% harder.
TEST(friction_is_the_same_along_a_diagonal) {
    for (Vec3 heading : {Vec3{1, 0, 0}, Vec3{0, 0, 1}, Vec3{1, 0, 1}, Vec3{-2, 0, 1}}) {
        World w;
        w.add(floor_box(0.4f));
        Body slab = Body::solid_box({0.5f, 0.1f, 0.5f}, 1, {0, 0.1f, 0});  // low and wide: it will not tip
        slab.friction = 0.4f;
        slab.restitution = 0;
        const Vec3 dir = heading.normalized();
        slab.vel = dir * 6;
        w.add(slab);
        run(w, 60);  // half a second
        const Vec3 v = w.bodies[1].vel;
        CHECK_NEAR(v.length(), 6 - 0.4 * kG * 0.5, 0.03);
        CHECK(dot(v.normalized(), dir) > 0.9999f);  // it does not curve
    }
}

// A solid sphere (I = 2/5 m r^2) sliding without spin is brought to rolling by friction at
// v = v0 / (1 + 2/5) = 5/7 v0, spinning about the horizontal axis perpendicular to its motion.
TEST(a_sliding_sphere_settles_into_rolling_at_five_sevenths) {
    World w;
    w.add(floor_box(0.5f));
    Body b = ball(0.5f, {0, 0.5f, 0});
    b.friction = 0.5f;
    b.vel = {7, 0, 0};
    w.add(b);
    run(w, 240);
    CHECK_VEC(w.bodies[1].vel, Vec3(5, 0, 0), 0.02);
    CHECK_VEC(w.bodies[1].w, Vec3(0, 0, -10), 0.05);  // moving +x on a floor below: clockwise seen from +z

    // The same along another heading: the spin axis turns with it.
    World w2;
    w2.add(floor_box(0.5f));
    Body c = ball(0.5f, {0, 0.5f, 0});
    c.friction = 0.5f;
    c.vel = {0, 0, 7};
    w2.add(c);
    run(w2, 240);
    CHECK_VEC(w2.bodies[1].vel, Vec3(0, 0, 5), 0.02);
    CHECK_VEC(w2.bodies[1].w, Vec3(10, 0, 0), 0.05);
}

// Rolling without slipping down a slope, a solid sphere accelerates at 5/7 g sin(theta).
TEST(a_sphere_rolls_down_a_slope_at_five_sevenths_g_sin_theta) {
    const Real theta = 0.3f;
    World w;
    Body ramp = Body::fixed_box({40, 0.5f, 20}, {0, 0, 0}, Quat::from_axis_angle({0, 0, 1}, theta));
    ramp.friction = 1;
    ramp.restitution = 0;
    w.add(ramp);
    Body b = ball(0.5f, Vec3{-std::sin(theta), std::cos(theta), 0} * 1.0f);
    b.friction = 1;
    w.add(b);
    const Vec3 down{-std::cos(theta), -std::sin(theta), 0};
    run(w, 120);
    const double v1 = double(dot(w.bodies[1].vel, down));
    run(w, 120);
    const double v2 = double(dot(w.bodies[1].vel, down));
    const double expected = 5.0 / 7.0 * kG * std::sin(double(theta));
    CHECK_NEAR(v2 - v1, expected, 0.02 * expected);
}

// ---- stacking and warm starting ------------------------------------------------------------------

namespace {

World tower(int crates, int sweeps, bool warm, int steps) {
    World w;
    w.solver.iterations = sweeps;
    w.solver.warm_starting = warm;
    w.add(floor_box());
    for (int i = 0; i < crates; ++i) w.add(crate({0, 0.5f + static_cast<Real>(i) * 1.001f, 0}));
    run(w, steps);
    return w;
}

}  // namespace

TEST(a_tower_of_crates_stands_in_3d) {
    World w = tower(4, 10, true, 600);
    for (int i = 1; i <= 4; ++i) {
        const Body& b = w.bodies[static_cast<size_t>(i)];
        CHECK_NEAR(b.pos.y, 0.5 + (i - 1), 0.02);
        CHECK_NEAR(b.pos.x, 0, 0.01);
        CHECK_NEAR(b.pos.z, 0, 0.01);
        CHECK_NEAR(double(b.vel.length()) + double(b.w.length()), 0, 1e-3);
    }
}

TEST(warm_starting_holds_stacks_that_the_cold_solver_drops) {
    CHECK(tower(4, 1, true, 600).bodies[4].pos.y > 3.4f);   // one sweep per step is enough when warm
    CHECK(tower(4, 1, false, 600).bodies[4].pos.y < 3.0f);  // and not when cold
    CHECK(tower(10, 4, true, 960).bodies[10].pos.y > 9.3f);
    CHECK(tower(10, 4, false, 960).bodies[10].pos.y < 8.0f);
}

TEST(a_pyramid_of_thirty_crates_stands) {
    World w;
    w.add(floor_box());
    for (int layer = 0; layer < 4; ++layer) {
        const int m = 4 - layer;  // 4x4, 3x3, 2x2, 1
        for (int i = 0; i < m; ++i)
            for (int k = 0; k < m; ++k)
                w.add(crate({(static_cast<Real>(i) - static_cast<Real>(m - 1) * 0.5f) * 1.02f, 0.5f + static_cast<Real>(layer) * 1.002f,
                             (static_cast<Real>(k) - static_cast<Real>(m - 1) * 0.5f) * 1.02f}));
    }
    CHECK(w.bodies.size() == 31);
    const auto start = positions(w);
    run(w, 8 * 120);
    const Settle s = settle(w, start);
    CHECK(s.moved < 0.03);
    CHECK(s.speed < 0.01);
}

// ---- solver plumbing -----------------------------------------------------------------------------

TEST(tangent_basis_is_orthonormal_and_right_handed_for_any_normal) {
    std::uint32_t seed = 3;
    auto rnd = [&] {
        seed = seed * 1664525u + 1013904223u;
        return static_cast<Real>(seed >> 8) / static_cast<Real>(1u << 24) * 2 - 1;
    };
    for (int i = 0; i < 500; ++i) {
        Vec3 n{rnd(), rnd(), rnd()};
        if (i < 6) {  // the axes themselves, both ways
            n = {};
            n[i / 2] = i % 2 ? -1.0f : 1.0f;
        }
        if (n.length() < 0.1f) continue;
        n = n.normalized();
        Vec3 t1, t2;
        tangent_basis(n, t1, t2);
        CHECK_NEAR(t1.length(), 1, 1e-4);
        CHECK_NEAR(t2.length(), 1, 1e-4);
        CHECK_NEAR(dot(t1, n), 0, 1e-4);
        CHECK_NEAR(dot(t2, n), 0, 1e-4);
        CHECK_NEAR(dot(t1, t2), 0, 1e-4);
        CHECK_VEC(cross(t1, t2), n, 1e-4);
    }
}

TEST(impulses_follow_the_same_pair_and_feature_id_in_3d) {
    auto pair_of = [](int a, int b, std::initializer_list<std::uint32_t> ids, Real n, Real t0, Real t1) {
        ContactPair p;
        p.a = a;
        p.b = b;
        for (std::uint32_t id : ids) {
            ContactPoint& cp = p.manifold.points[p.manifold.count++];
            cp.id = id;
            cp.normal_impulse = n;
            cp.tangent_impulse[0] = t0;
            cp.tangent_impulse[1] = t1;
        }
        return p;
    };
    std::vector<ContactPair> previous = {pair_of(1, 2, {5, 7}, 3.0f, 0.5f, -0.25f), pair_of(1, 3, {5}, 9.0f, 0, 0)};
    previous[0].manifold.points[1].normal_impulse = 4.0f;
    std::vector<ContactPair> current = {pair_of(1, 2, {7, 9}, 99, 99, 99), pair_of(1, 3, {6}, 99, 99, 99), pair_of(2, 3, {5}, 0, 0, 0)};
    transfer_impulses(previous, current);
    CHECK_NEAR(current[0].manifold.points[0].normal_impulse, 4.0, 1e-6);  // id 7 inherited, all three components
    CHECK_NEAR(current[0].manifold.points[0].tangent_impulse[0], 0.5, 1e-6);
    CHECK_NEAR(current[0].manifold.points[0].tangent_impulse[1], -0.25, 1e-6);
    CHECK_NEAR(current[0].manifold.points[1].normal_impulse, 0, 1e-6);    // id 9 is new: stale values cleared
    CHECK_NEAR(current[0].manifold.points[1].tangent_impulse[1], 0, 1e-6);
    CHECK_NEAR(current[1].manifold.points[0].normal_impulse, 0, 1e-6);    // different feature
    CHECK_NEAR(current[2].manifold.points[0].normal_impulse, 0, 1e-6);    // different pair
}

TEST(the_same_scene_steps_identically_twice) {
    auto build = [] {
        World w;
        w.add(floor_box());
        for (int i = 0; i < 3; ++i) w.add(crate({0.1f * static_cast<Real>(i), 0.5f + static_cast<Real>(i) * 1.2f, 0.05f * static_cast<Real>(i)}));
        w.add(ball(0.4f, {0.3f, 6, 0.2f}, 0.5f));
        return w;
    };
    World a = build(), b = build();
    run(a, 400);
    run(b, 400);
    for (size_t i = 0; i < a.bodies.size(); ++i) {
        CHECK(a.bodies[i].pos.x == b.bodies[i].pos.x && a.bodies[i].pos.y == b.bodies[i].pos.y && a.bodies[i].pos.z == b.bodies[i].pos.z);
        CHECK(a.bodies[i].q.w == b.bodies[i].q.w && a.bodies[i].w.x == b.bodies[i].w.x);
    }
}

TEST(truncate_forgets_contacts_in_3d) {
    World w;
    w.add(floor_box());
    w.add(crate({0, 0.5f, 0}));
    run(w, 60);
    CHECK(!w.contacts().empty());
    w.truncate(1);
    CHECK(w.bodies.size() == 1 && w.contacts().empty());
    w.add(crate({0, 3, 0}));
    w.step(kDt);
    CHECK(w.contacts().empty());
}

// ---- gaps found by breaking the solver on purpose ------------------------------------------------

// A contact impulse is equal and opposite at one point, so whatever the shapes and however they are
// turned, one solve must leave the pair's total momentum and total angular momentum (about the origin)
// unchanged. This is the test that the solver uses each body's inertia in the WORLD frame: with tilted,
// brick-shaped bodies the body-frame tensor is a different matrix.
TEST(a_solve_conserves_momentum_and_angular_momentum_for_tilted_bricks) {
    std::uint32_t seed = 17;
    auto rnd = [&](Real lo, Real hi) {
        seed = seed * 1664525u + 1013904223u;
        return lo + (hi - lo) * static_cast<Real>(seed >> 8) / static_cast<Real>(1u << 24);
    };
    auto rvec = [&](Real e) { return Vec3{rnd(-e, e), rnd(-e, e), rnd(-e, e)}; };
    int solved = 0;
    for (int trial = 0; trial < 600; ++trial) {
        std::vector<Body> bodies;
        for (int k = 0; k < 2; ++k) {
            Body b = Body::solid_box({rnd(0.3f, 1.0f), rnd(0.2f, 0.6f), rnd(0.1f, 0.4f)}, rnd(0.5f, 3), rvec(0.6f),
                                     Quat::from_axis_angle(rvec(1), rnd(0, 3)));
            b.vel = rvec(2);
            b.w = rvec(3);
            b.restitution = 0.5f;
            bodies.push_back(b);
        }
        ContactPair pair;
        pair.a = 0;
        pair.b = 1;
        if (!collide(bodies[0], bodies[1], pair.manifold)) continue;
        std::vector<ContactPair> contacts = {pair};
        auto momentum = [&] { return bodies[0].vel * bodies[0].mass + bodies[1].vel * bodies[1].mass; };
        auto angular = [&] {
            Vec3 l;
            for (const Body& b : bodies) l += cross(b.pos, b.vel) * b.mass + b.angular_momentum();
            return l;
        };
        const Vec3 p0 = momentum(), l0 = angular();
        const Real scale = 1 + p0.length() + l0.length();
        solve_contacts(bodies, contacts, kDt, SolverSettings{});
        CHECK_VEC(momentum(), p0, 2e-4f * scale);
        CHECK_VEC(angular(), l0, 5e-4f * scale);
        ++solved;
    }
    CHECK(solved > 150);
}

// A crate held on a slope by static friction: the stored friction impulses, added up, are exactly what
// cancels gravity's pull along the slope each step, m g sin(theta) dt, pointing uphill. Tilting the slope
// about x instead of z puts that pull along the other tangent axis, so both stored components are checked.
TEST(stored_friction_impulses_hold_a_crate_on_a_slope) {
    for (int about = 0; about < 2; ++about) {
        const Vec3 axis = about == 0 ? Vec3{0, 0, 1} : Vec3{1, 0, 0};
        const Quat tilt = Quat::from_axis_angle(axis, kTheta);
        World w;
        Body ramp = Body::fixed_box({20, 0.5f, 20}, {0, 0, 0}, tilt);
        ramp.friction = 0.6f;
        ramp.restitution = 0;
        w.add(ramp);
        const Vec3 normal = rotate(tilt, {0, 1, 0});
        Body slider = Body::solid_box({0.5f, 0.5f, 0.5f}, 1, normal * 1.0f, tilt);
        slider.friction = 0.6f;
        slider.restitution = 0;
        w.add(slider);
        run(w, 360);
        CHECK_NEAR(w.bodies[1].vel.length(), 0, 0.01);
        CHECK(w.contacts().size() == 1);
        if (w.contacts().empty()) continue;  // already failed; do not index into nothing
        const Manifold& m = w.contacts()[0].manifold;
        Vec3 t1, t2, friction;
        tangent_basis(m.normal, t1, t2);
        for (int k = 0; k < m.count; ++k) friction += t1 * m.points[k].tangent_impulse[0] + t2 * m.points[k].tangent_impulse[1];
        // The force on the crate (body b) is +friction. Gravity along the slope is g minus its normal part.
        const Vec3 g{0, static_cast<Real>(-kG), 0};
        const Vec3 along = g - normal * dot(g, normal);
        CHECK_VEC(friction, along * -kDt, 0.03 * kG * std::sin(double(kTheta)) * double(kDt));
        CHECK(friction.length() > 0.01f);
    }
}

TEST(a_body_that_starts_inside_the_floor_is_pushed_out) {
    World w;
    w.add(floor_box());
    w.add(crate({0, 0.44f, 0}));  // 6 cm into the floor
    run(w, 360);
    CHECK(w.bodies[1].pos.y > 0.49f && w.bodies[1].pos.y < 0.502f);  // back to within the slop of resting height
    CHECK_NEAR(w.bodies[1].vel.length(), 0, 0.01);
}

// Friction between two materials is the geometric mean, so it does not matter which body is which.
TEST(friction_mixes_both_materials_whichever_is_first) {
    auto speed_after = [](Real floor_mu, Real slab_mu) {
        World w;
        w.add(floor_box(floor_mu));
        Body slab = Body::solid_box({0.5f, 0.1f, 0.5f}, 1, {0, 0.1f, 0});
        slab.friction = slab_mu;
        slab.restitution = 0;
        slab.vel = {6, 0, 0};
        w.add(slab);
        run(w, 60);
        return double(w.bodies[1].vel.x);
    };
    const double expected = 6 - std::sqrt(0.9 * 0.1) * kG * 0.5;  // mu = sqrt(0.9 x 0.1) = 0.3
    CHECK_NEAR(speed_after(0.9f, 0.1f), expected, 0.03);
    CHECK_NEAR(speed_after(0.1f, 0.9f), expected, 0.03);
}

// The World passes its gyroscopic setting on: a brick tumbling in empty space keeps its angular momentum
// pointing the same way by default, and visibly does not when the term is switched off.
TEST(the_world_integrates_tumbling_bodies_with_its_gyroscopic_setting) {
    auto momentum_swing = [](Gyroscopic mode) {
        World w;
        w.gravity = {};
        w.gyroscopic = mode;
        Body b = Body::solid_box({1.0f, 0.6f, 0.2f}, 1, {0, 0, 0});
        b.w = {3, 2, 1};
        w.add(b);
        const Vec3 l0 = w.bodies[0].angular_momentum();
        double worst = 0;
        for (int i = 0; i < 5 * 240; ++i) {
            w.step(1.0f / 240.0f);
            const Vec3 l = w.bodies[0].angular_momentum();
            worst = std::max(worst, std::acos(std::min(1.0, double(dot(l, l0)) / (double(l.length()) * double(l0.length())))));
        }
        return worst;
    };
    CHECK(momentum_swing(Gyroscopic::Implicit) < 0.03);
    CHECK(momentum_swing(Gyroscopic::Off) > 0.3);
    World defaults;
    CHECK(defaults.gyroscopic == Gyroscopic::Implicit);
}

// Friction is warm started like the normal impulse. With a single sweep per step a crate on a slope stays
// dead still because each step begins from the friction that held it last step; without that it creeps
// downhill at 4 to 10 mm/s (measured), never quite coming to rest.
TEST(friction_is_warm_started_too) {
    for (int about = 0; about < 2; ++about) {
        const Quat tilt = Quat::from_axis_angle(about == 0 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}, kTheta);
        World w;
        w.solver.iterations = 1;
        Body ramp = Body::fixed_box({20, 0.5f, 20}, {0, 0, 0}, tilt);
        ramp.friction = 0.6f;
        ramp.restitution = 0;
        w.add(ramp);
        Body slider = Body::solid_box({0.5f, 0.5f, 0.5f}, 1, rotate(tilt, {0, 1, 0}) * 1.0f, tilt);
        slider.friction = 0.6f;
        slider.restitution = 0;
        w.add(slider);
        const Vec3 start = w.bodies[1].pos;
        run(w, 480);
        CHECK(distance(w.bodies[1].pos, start) < 0.006f);  // 4 mm of settling, then nothing
        CHECK_NEAR(w.bodies[1].vel.length(), 0, 0.001);
    }
}
