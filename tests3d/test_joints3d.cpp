#include "test.hpp"
#include <phys3d/world.hpp>

#include <phys/world.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace phys3d;

namespace {

constexpr Real kG = 9.81f;
constexpr Real kDt = 1.0f / 240.0f;

void run(World& w, int steps, Real dt = kDt) {
    for (int i = 0; i < steps; ++i) w.step(dt);
}
// Distance between a joint's two anchors: how far a pin joint has come apart.
Real gap(const World& w, int joint) {
    const Joint& j = w.joints[static_cast<size_t>(joint)];
    return distance(j.world_anchor_a(w.bodies), j.world_anchor_b(w.bodies));
}
Real total_energy(const World& w) {
    Real e = 0;
    for (const Body& b : w.bodies)
        if (b.type == BodyType::Dynamic) e += b.kinetic_energy() - b.mass * dot(w.gravity, b.pos);
    return e;
}
// Two tumbling bodies in arbitrary poses, for checking rows against finite differences.
std::vector<Body> two_tumblers() {
    std::vector<Body> bodies;
    Body a = Body::solid_box({0.3f, 0.5f, 0.2f}, 2, {0.2f, 1.0f, -0.4f}, Quat::from_axis_angle({1, 2, 3}, 0.7f));
    a.vel = {0.3f, -0.6f, 0.9f};
    a.w = {1.1f, -0.7f, 0.5f};
    Body b = Body::solid_box({0.4f, 0.2f, 0.6f}, 1, {1.3f, 0.6f, 0.5f}, Quat::from_axis_angle({-2, 1, 1}, -1.1f));
    b.vel = {-0.8f, 0.2f, 0.4f};
    b.w = {-0.4f, 0.9f, 1.3f};
    bodies.push_back(a);
    bodies.push_back(b);
    return bodies;
}
// Each row claims "J v is how fast my error is changing". Check that against the error actually measured a
// moment before and a moment after, with the bodies coasting on their velocities.
void check_rows_against_finite_differences(const std::vector<Body>& bodies, Joint joint, size_t expected_rows) {
    const Real h = 1e-3f;
    auto rows_at = [&](Real t) {
        std::vector<Body> moved = bodies;
        for (Body& b : moved) b.integrate_position(t);
        JointSolver solver;
        solver.prepare(moved, std::span<Joint>(&joint, 1), kDt, SolverSettings{});
        std::vector<std::pair<Real, Real>> out;  // (error, J v)
        for (const JointRow& r : solver.rows()) out.emplace_back(r.error, r.velocity());
        return out;
    };
    const auto before = rows_at(-h), now = rows_at(0), after = rows_at(h);
    CHECK(now.size() == expected_rows);
    if (before.size() != expected_rows || now.size() != expected_rows || after.size() != expected_rows) return;
    for (size_t i = 0; i < now.size(); ++i) CHECK_NEAR((after[i].first - before[i].first) / (2 * h), now[i].second, 5e-3);
}

}  // namespace

// ---- the rows ------------------------------------------------------------------------------------

TEST(distance_row_matches_the_rate_of_change_of_the_distance) {
    const std::vector<Body> bodies = two_tumblers();
    check_rows_against_finite_differences(bodies, Joint::distance(bodies, 0, 1, {0.4f, 1.2f, -0.3f}, {1.1f, 0.5f, 0.8f}), 1);
    check_rows_against_finite_differences(bodies, Joint::distance(bodies, -1, 1, {0, 3, 0}, {1.1f, 0.5f, 0.8f}), 1);
}

TEST(ball_rows_match_the_rate_the_anchors_separate) {
    const std::vector<Body> bodies = two_tumblers();
    check_rows_against_finite_differences(bodies, Joint::ball(bodies, 0, 1, {0.7f, 0.9f, 0.1f}), 3);
    check_rows_against_finite_differences(bodies, Joint::ball(bodies, 1, -1, {0.7f, 0.9f, 0.1f}), 3);
}

TEST(hinge_rows_match_the_rate_the_axes_and_the_angle_change) {
    const std::vector<Body> bodies = two_tumblers();
    Joint hinge = Joint::hinge(bodies, 0, 1, {0.7f, 0.9f, 0.1f}, {1, -2, 0.5f});
    check_rows_against_finite_differences(bodies, hinge, 5);  // 3 for the point, 2 for the axis
    // Both limits violated, so both limit rows report an error that changes with the angle.
    hinge.enable_limit = true;
    hinge.lower = 1;
    hinge.upper = -1;
    check_rows_against_finite_differences(bodies, hinge, 7);
}

TEST(hinge_angle_reads_back_a_known_turn) {
    const Vec3 axis = Vec3{1, 2, 3}.normalized();
    for (Real turn : {-3.0f, -1.0f, 0.0f, 0.5f, 3.0f}) {
        std::vector<Body> bodies;
        bodies.push_back(Body::solid_box({0.2f, 0.3f, 0.4f}, 1, {1, 2, 3}, Quat::from_axis_angle({3, -1, 2}, 0.9f)));
        const Joint hinge = Joint::hinge(bodies, -1, 0, {1, 2, 3}, axis);
        CHECK_NEAR(hinge.hinge_angle(bodies), 0, 1e-5);
        bodies[0].q = Quat::from_axis_angle(axis, turn) * bodies[0].q;
        CHECK_NEAR(hinge.hinge_angle(bodies), turn, 1e-4);
        // Measured the other way round, with the body as `a`, the same turn reads negative.
        bodies[0].q = Quat::from_axis_angle({3, -1, 2}, 0.9f);
        const Joint flipped = Joint::hinge(bodies, 0, -1, {1, 2, 3}, axis);
        bodies[0].q = Quat::from_axis_angle(axis, turn) * bodies[0].q;
        CHECK_NEAR(flipped.hinge_angle(bodies), -turn, 1e-4);
    }
}

// ---- ball and socket -----------------------------------------------------------------------------

TEST(a_ball_jointed_pendulum_swings_with_the_textbook_period) {
    // A solid sphere hung from the world by a ball joint a length L from its centre: a compound pendulum,
    // T = 2 pi sqrt((I + m L^2) / (m g L)) with I = 2/5 m r^2.
    const Real L = 1, r = 0.1f, start = 0.1f;
    World w;
    const int bob = w.add(Body::solid_sphere(r, 1, {L * std::sin(start), -L * std::cos(start), 0}));
    const int pin = w.add_joint(Joint::ball(w.bodies, -1, bob, {0, 0, 0}));
    const Real e0 = total_energy(w);
    const double expected = 2 * kPi * std::sqrt((0.4 * r * r + L * L) / (kG * L));

    std::vector<double> crossings;  // times the bob passes the bottom moving in +x
    Real worst_gap = 0, worst_energy = 0, last_x = w.bodies[0].pos.x;
    for (int i = 1; i <= 240 * 10; ++i) {
        w.step(kDt);
        const Real x = w.bodies[0].pos.x;
        if (last_x < 0 && x >= 0) crossings.push_back((i - 1 + last_x / (last_x - x)) * static_cast<double>(kDt));
        last_x = x;
        worst_gap = std::max(worst_gap, gap(w, pin));
        worst_energy = std::max(worst_energy, std::fabs(total_energy(w) - e0));
    }
    CHECK(crossings.size() >= 4);
    if (crossings.size() < 4) return;
    CHECK_NEAR((crossings[3] - crossings[0]) / 3, expected, expected * 0.005);
    CHECK(worst_gap < 1e-4f);
    // The swing holds 0.5 m g L start^2 of energy; it must keep nearly all of it over five swings.
    CHECK(worst_energy < 0.03f * 0.5f * w.bodies[0].mass * kG * L * start * start);
    CHECK(std::fabs(w.bodies[0].pos.z) < 1e-4f);  // and stays in its plane
}

TEST(a_ball_joint_leaves_all_three_rotations_free) {
    // No gravity, pinned at its own centre: the joint must not touch the spin at all. (A sphere, so there
    // is no gyroscopic wobble to confuse the picture.)
    World w;
    w.gravity = {};
    const int ball = w.add(Body::solid_sphere(0.3f, 1, {1, 2, 3}));
    w.add_joint(Joint::ball(w.bodies, -1, ball, {1, 2, 3}));
    w.bodies[0].w = {1.5f, -2.5f, 0.7f};
    run(w, 480);
    CHECK_NEAR(w.bodies[0].w.x, 1.5, 1e-4);
    CHECK_NEAR(w.bodies[0].w.y, -2.5, 1e-4);
    CHECK_NEAR(w.bodies[0].w.z, 0.7, 1e-4);
    CHECK(distance(w.bodies[0].pos, Vec3{1, 2, 3}) < 1e-5f);
}

TEST(a_pendulum_swung_in_a_circle_keeps_its_angular_momentum_about_the_vertical) {
    // Gravity and the pin's force both pass through, or run parallel to, the vertical line through the
    // pivot, so neither can change the angular momentum about it. This is the test a planar engine cannot
    // run. (A fine timestep: see the next test for what a coarse one costs.)
    World w;
    const Real dt = 1.0f / 960.0f;
    const int bob = w.add(Body::solid_sphere(0.1f, 1, {0.6f, -0.8f, 0}));
    w.add_joint(Joint::ball(w.bodies, -1, bob, {0, 0, 0}));
    w.bodies[0].vel = {0, 0, 1.2f};  // sideways: it will go round, rising and falling
    auto about_vertical = [&] { return (cross(w.bodies[0].pos, w.bodies[0].vel) * w.bodies[0].mass + w.bodies[0].angular_momentum()).y; };
    const Real l0 = about_vertical(), e0 = total_energy(w), kinetic0 = w.bodies[0].kinetic_energy();
    Real lowest = 0, highest = -10, worst_l = 0, worst_radius = 0;
    for (int i = 0; i < 960 * 3; ++i) {
        w.step(dt);
        lowest = std::min(lowest, w.bodies[0].pos.y);
        highest = std::max(highest, w.bodies[0].pos.y);
        worst_l = std::max(worst_l, std::fabs(about_vertical() - l0));
        worst_radius = std::max(worst_radius, std::fabs(w.bodies[0].pos.length() - 1));
    }
    CHECK(worst_l < 0.01f * std::fabs(l0));
    CHECK_NEAR(total_energy(w), e0, 0.03 * kinetic0);  // it loses about 2%
    CHECK(highest - lowest > 0.1f);  // it really did rise and fall
    CHECK(worst_radius < 1e-4f);     // on the sphere of radius 1 about the pivot
}

TEST(spin_about_an_off_centre_pin_is_lost_in_proportion_to_the_timestep) {
    // Not a property to be proud of, but one to know. Each step the body moves along the tangent, the
    // position pass pulls it back onto the circle, and the velocity is then projected onto the new tangent:
    // a little speed is lost every step, about (w dt)^2 of the energy. It is the price of a first-order
    // integrator, the 2D engine pays the same, and it shrinks with the step.
    auto speed_left = [](Real hz) {
        World w;
        w.gravity = {};
        const int door = w.add(Body::solid_box({0.5f, 0.3f, 0.05f}, 1, {0.5f, 0, 0}));
        w.add_joint(Joint::hinge(w.bodies, -1, door, {0, 0, 0}, {0, 1, 0}));
        w.bodies[0].w = {0, 3, 0};
        w.bodies[0].vel = cross(Vec3{0, 3, 0}, Vec3{0.5f, 0, 0});  // already circling the pin
        run(w, static_cast<int>(hz) * 5, 1 / hz);
        return w.bodies[0].w.y;
    };
    const Real coarse = 3 - speed_left(240), fine = 3 - speed_left(960);
    CHECK(coarse > 0 && coarse < 0.25f);  // 6% in five seconds at 240 Hz
    CHECK(fine > 0 && fine < 0.3f * coarse);
}

TEST(two_free_bodies_joined_by_a_joint_keep_their_total_momentum) {
    // A joint is an internal force: whatever it does to one body it does the opposite to the other, at the
    // same point. Linear and angular momentum of the pair (about the origin) cannot change.
    for (int kind = 0; kind < 3; ++kind) {
        World w;
        w.gravity = {};
        for (const Body& b : two_tumblers()) w.add(b);
        const Vec3 anchor{0.7f, 0.9f, 0.1f};
        w.add_joint(kind == 0   ? Joint::distance(w.bodies, 0, 1, {0.4f, 1.2f, -0.3f}, {1.1f, 0.5f, 0.8f})
                    : kind == 1 ? Joint::ball(w.bodies, 0, 1, anchor)
                                : Joint::hinge(w.bodies, 0, 1, anchor, {1, -2, 0.5f}));
        auto momentum = [&](Vec3& linear, Vec3& angular) {
            linear = angular = {};
            for (const Body& b : w.bodies) {
                linear += b.vel * b.mass;
                angular += cross(b.pos, b.vel) * b.mass + b.angular_momentum();
            }
        };
        Vec3 p0, l0, p1, l1;
        momentum(p0, l0);
        const Real e0 = total_energy(w);
        run(w, 480);
        momentum(p1, l1);
        CHECK(distance(p0, p1) < 1e-3f * p0.length());
        CHECK(distance(l0, l1) < 0.02f * l0.length());  // the gyroscopic step is not exactly momentum-conserving
        CHECK(total_energy(w) <= e0 * 1.001f);          // and a joint never adds energy
        CHECK(total_energy(w) > e0 * 0.5f);             // nor quietly eats it all
    }
}

TEST(a_resting_pendulum_stores_exactly_the_impulse_that_carries_its_weight) {
    World w;
    const int bob = w.add(Body::solid_sphere(0.2f, 3, {0, -1, 0}));
    const int pin = w.add_joint(Joint::ball(w.bodies, -1, bob, {0, 0, 0}));
    run(w, 240);
    const Joint& j = w.joints[static_cast<size_t>(pin)];
    const Real weight_impulse = w.bodies[0].mass * kG * kDt;
    CHECK_NEAR(j.impulses[0], 0, 1e-4 * weight_impulse);
    CHECK_NEAR(j.impulses[1], weight_impulse, 1e-3 * weight_impulse);  // +y on the bob, each step
    CHECK_NEAR(j.impulses[2], 0, 1e-4 * weight_impulse);
}

TEST(a_chain_carrying_a_heavy_weight_stays_together) {
    // The configuration that tore 2D chains apart when position error was fed back as velocity.
    World w;
    const int links = 8;
    int previous = -1;
    for (int i = 0; i < links; ++i) {
        const int link = w.add(Body::solid_box({0.2f, 0.04f, 0.04f}, 1, {0.2f + 0.4f * static_cast<Real>(i), 0, 0}));
        w.add_joint(Joint::ball(w.bodies, previous, link, {0.4f * static_cast<Real>(i), 0, 0}));
        previous = link;
    }
    const int weight = w.add(Body::solid_sphere(0.25f, 1, {0.4f * links + 0.25f, 0, 0}));  // 26x a link
    w.add_joint(Joint::ball(w.bodies, previous, weight, {0.4f * links, 0, 0}));
    const Real e0 = total_energy(w);
    Real worst_gap = 0, worst_energy = e0;
    for (int i = 0; i < 240 * 6; ++i) {
        w.step(kDt);
        for (int j = 0; j <= links; ++j) worst_gap = std::max(worst_gap, gap(w, j));
        worst_energy = std::max(worst_energy, total_energy(w));
    }
    CHECK(worst_gap < 0.02f);  // under 1% of the chain's length, at the bottom of a 3 m swing
    CHECK(worst_energy <= e0 + 1e-4f);  // it started at rest: it may only lose energy
    CHECK(w.bodies[static_cast<size_t>(weight)].pos.y < -1);   // and it did swing down
}

TEST(the_position_pass_closes_a_gap_without_creating_motion) {
    World w;
    w.gravity = {};
    const int a = w.add(Body::solid_box({0.2f, 0.2f, 0.2f}, 1, {0, 0, 0}));
    const int b = w.add(Body::solid_box({0.2f, 0.2f, 0.2f}, 1, {0.4f, 0, 0}));
    const int pin = w.add_joint(Joint::hinge(w.bodies, a, b, {0.2f, 0, 0}, {0, 0, 1}));
    w.bodies[1].pos += Vec3{0.05f, 0.08f, -0.03f};  // pulled apart by hand...
    w.bodies[1].q = Quat::from_axis_angle({1, 0, 0}, 0.1f);  // ...and twisted off the hinge axis
    run(w, 20);
    CHECK(gap(w, pin) < 1e-4f);
    const Vec3 axis_b = rotate(w.bodies[1].q, w.joints[0].local_axis_b);
    CHECK(cross(w.joints[0].world_axis(w.bodies), axis_b).length() < 1e-3f);
    for (const Body& body : w.bodies) CHECK(body.vel.length() == 0 && body.w.length() == 0);
}

// ---- distance ------------------------------------------------------------------------------------

TEST(a_rigid_distance_joint_holds_its_length) {
    World w;
    const int bob = w.add(Body::solid_box({0.1f, 0.2f, 0.15f}, 1, {1, 0, 0}));
    const int rod = w.add_joint(Joint::distance(w.bodies, -1, bob, {0, 0, 0}, {0.9f, 0.2f, 0.15f}));  // a corner
    const Real length = w.joints[0].length;
    w.bodies[0].vel = {0, 0, 2};
    Real worst = 0;
    for (int i = 0; i < 240 * 4; ++i) {
        w.step(kDt);
        worst = std::max(worst, std::fabs(gap(w, rod) - length));
    }
    CHECK(worst < 2e-3f);
}

TEST(a_soft_distance_joint_is_a_spring_of_the_asked_frequency) {
    World w;
    w.gravity = {};
    const int bob = w.add(Body::solid_sphere(0.1f, 1, {1, 0, 0}));
    Joint spring = Joint::distance(w.bodies, -1, bob, {0, 0, 0}, {1, 0, 0});
    spring.frequency_hz = 2;
    spring.damping_ratio = 0;
    w.add_joint(spring);
    w.bodies[0].vel = {0.5f, 0, 0};
    std::vector<double> crossings;
    Real last = 0;
    for (int i = 1; i <= 240 * 3; ++i) {
        w.step(kDt);
        const Real x = w.bodies[0].pos.x - 1;
        if (last < 0 && x >= 0) crossings.push_back((i - 1 + last / (last - x)) * static_cast<double>(kDt));
        last = x;
    }
    CHECK(crossings.size() >= 4);
    if (crossings.size() < 4) return;
    CHECK_NEAR((crossings[3] - crossings[0]) / 3, 0.5, 0.01);
}

TEST(a_damped_spring_settles_and_an_undamped_one_does_not) {
    auto speed_after = [](Real damping_ratio) {
        World w;
        w.gravity = {};
        w.allow_sleep = false;
        const int bob = w.add(Body::solid_sphere(0.1f, 1, {1, 0, 0}));
        Joint spring = Joint::distance(w.bodies, -1, bob, {0, 0, 0}, {1, 0, 0});
        spring.frequency_hz = 2;
        spring.damping_ratio = damping_ratio;
        w.add_joint(spring);
        w.bodies[0].vel = {0.5f, 0, 0};
        Real peak = 0;
        for (int i = 0; i < 240 * 3; ++i) {
            w.step(kDt);
            if (i >= 240 * 2) peak = std::max(peak, w.bodies[0].vel.length());
        }
        return peak;
    };
    CHECK(speed_after(1) < 1e-3f);
    // (The implicit spring is not perfectly lossless: at 240 Hz about half the speed is left by now.)
    CHECK(speed_after(0) > 0.2f);
}

TEST(a_lightly_damped_spring_decays_at_the_textbook_rate) {
    // Damping ratio z: each swing's peak is exp(-2 pi z / sqrt(1 - z^2)) of the one before.
    World w;
    w.gravity = {};
    w.allow_sleep = false;
    const Real dt = 1.0f / 960.0f, zeta = 0.2f;
    const int bob = w.add(Body::solid_sphere(0.1f, 1, {1, 0, 0}));
    Joint spring = Joint::distance(w.bodies, -1, bob, {0, 0, 0}, {1, 0, 0});
    spring.frequency_hz = 2;
    spring.damping_ratio = zeta;
    w.add_joint(spring);
    w.bodies[0].vel = {0.5f, 0, 0};
    std::vector<Real> peaks;  // successive maxima of the stretch
    Real before = 0, last = 0;
    for (int i = 0; i < 960 * 2; ++i) {
        w.step(dt);
        const Real x = w.bodies[0].pos.x - 1;
        if (last > before && last > x && last > 0) peaks.push_back(last);
        before = last;
        last = x;
    }
    CHECK(peaks.size() >= 3);
    if (peaks.size() < 3) return;
    const double expected = std::exp(-2 * kPi * zeta / std::sqrt(1 - zeta * zeta));  // 0.277
    CHECK_NEAR(peaks[1] / peaks[0], expected, 0.02);
    CHECK_NEAR(peaks[2] / peaks[1], expected, 0.02);
}

TEST(a_badly_violated_limit_is_eased_back_not_snapped) {
    // One correction turns a body at most 0.14 rad (8 degrees), so a joint that starts a long way out does
    // not teleport.
    World w;
    w.gravity = {};
    w.solver.joint_position_iterations = 1;
    const int wheel = w.add(Body::solid_box({0.3f, 0.3f, 0.1f}, 1, {0, 0, 0}));
    Joint hinge = Joint::hinge(w.bodies, -1, wheel, {0, 0, 0}, {0, 0, 1});
    hinge.enable_limit = true;
    hinge.lower = -0.1f;
    hinge.upper = 0.1f;
    w.add_joint(hinge);
    w.bodies[0].q = Quat::from_axis_angle({0, 0, 1}, 1.1f);  // a whole radian past the upper limit
    w.step(kDt);
    CHECK_NEAR(w.joints[0].hinge_angle(w.bodies), 1.1 - 0.14, 2e-3);
    CHECK(w.bodies[0].w.length() == 0);
    run(w, 20);
    CHECK_NEAR(w.joints[0].hinge_angle(w.bodies), 0.1, 1e-3);  // and it does get there
}

// ---- hinge ---------------------------------------------------------------------------------------

TEST(a_hinge_allows_turning_about_its_axis_and_nothing_else) {
    World w;
    w.gravity = {};
    const Vec3 axis = Vec3{1, 1, 0}.normalized();
    const int door = w.add(Body::solid_box({0.5f, 0.3f, 0.05f}, 1, {0.5f, 0, 0}));
    const int pin = w.add_joint(Joint::hinge(w.bodies, -1, door, {0, 0, 0}, axis));
    // Knock it every which way.
    w.bodies[0].vel = {1, -2, 3};
    w.bodies[0].w = Vec3{4, -3, 5};
    Real worst_gap = 0, worst_tilt = 0;
    for (int i = 0; i < 240 * 4; ++i) {
        w.step(kDt);
        worst_gap = std::max(worst_gap, gap(w, pin));
        worst_tilt = std::max(worst_tilt, cross(axis, rotate(w.bodies[0].q, w.joints[0].local_axis_b)).length());
    }
    CHECK(worst_gap < 2e-3f);
    CHECK(worst_tilt < 5e-3f);  // under a third of a degree
    const Vec3 spin = w.bodies[0].w;
    CHECK(cross(spin, axis).length() < 1e-3f * spin.length());  // what is left is spin about the axis
    CHECK(spin.length() > 0.5f);                                 // and there is some
}

TEST(a_free_hinge_does_not_slow_the_spin_about_its_axis) {
    // Pinned through its own centre, so nothing has to be pulled round a circle.
    World w;
    w.gravity = {};
    const int wheel = w.add(Body::solid_box({0.5f, 0.05f, 0.3f}, 1, {1, 2, 3}));
    w.add_joint(Joint::hinge(w.bodies, -1, wheel, {1, 2, 3}, {0, 1, 0}));
    w.bodies[0].w = {0, 3, 0};
    run(w, 240 * 5);
    CHECK_NEAR(w.bodies[0].w.y, 3, 1e-4);
}

TEST(a_hinged_arm_matches_the_2d_engine) {
    // The same swinging arm, pinned at one end, in both engines: same mass, same inertia about the hinge
    // axis. Turning about a principal axis has no gyroscopic term, so the two must agree.
    const Real mass = 2, half_length = 0.6f, inertia = mass * (half_length * half_length + 0.1f * 0.1f) / 3;
    World w3;
    Body arm3 = Body::dynamic(mass, Mat3::diagonal(inertia * 0.5f, inertia * 0.7f, inertia), {half_length, 0, 0});
    w3.add(arm3);
    w3.add_joint(Joint::hinge(w3.bodies, -1, 0, {0, 0, 0}, {0, 0, 1}));

    phys::World w2;
    phys::Body arm2(phys::Shape::make_polygon(phys::Polygon::box(half_length, 0.1f)), {half_length, 0}, 0);
    arm2.mass = mass;
    arm2.inv_mass = 1 / mass;
    arm2.inertia = inertia;
    arm2.inv_inertia = 1 / inertia;
    w2.add(arm2);
    w2.add_joint(phys::Joint::revolute(w2.bodies, -1, 0, {0, 0}));

    for (int i = 0; i < 240 * 3; ++i) {
        w3.step(kDt);
        w2.step(kDt);
    }
    CHECK(std::fabs(w3.bodies[0].pos.y) > 0.05f || std::fabs(w3.bodies[0].pos.x) < half_length - 0.05f);  // it moved
    CHECK_NEAR(w3.bodies[0].pos.x, w2.bodies[0].pos.x, 2e-3);
    CHECK_NEAR(w3.bodies[0].pos.y, w2.bodies[0].pos.y, 2e-3);
    CHECK_NEAR(w3.bodies[0].w.z, w2.bodies[0].w, 5e-3);
    CHECK(std::fabs(w3.bodies[0].pos.z) < 1e-5f);
}

TEST(a_torque_limited_motor_spins_a_wheel_up_at_torque_over_inertia) {
    // From rest with no gravity, a motor at its torque limit gives  w = torque * t / I.
    World w;
    w.gravity = {};
    const Vec3 half{0.4f, 0.3f, 0.1f};
    const int wheel = w.add(Body::solid_box(half, 1, {0, 0, 0}));
    Joint hinge = Joint::hinge(w.bodies, -1, wheel, {0, 0, 0}, {0, 0, 1});
    hinge.enable_motor = true;
    hinge.motor_speed = 1000;  // never reached
    hinge.max_motor = 0.02f;
    w.add_joint(hinge);
    const Real about_axis = w.bodies[0].mass * (half.x * half.x + half.y * half.y) / 3;
    run(w, 240);
    CHECK_NEAR(w.bodies[0].w.z, 0.02 / about_axis, 0.002 * 0.02 / about_axis);
    CHECK(w.bodies[0].w.z < 100);  // sanity: well short of the target speed
}

TEST(a_strong_motor_holds_its_speed_and_reverses) {
    World w;
    w.gravity = {};
    const int a = w.add(Body::solid_box({0.3f, 0.3f, 0.1f}, 1, {0, 0, 0}));
    const int b = w.add(Body::solid_box({0.3f, 0.3f, 0.1f}, 1, {0, 0, 0.2f}));
    Joint hinge = Joint::hinge(w.bodies, a, b, {0, 0, 0.1f}, {0, 0, 1});
    hinge.enable_motor = true;
    hinge.motor_speed = 2;
    hinge.max_motor = 100;
    const int motor = w.add_joint(hinge);
    run(w, 60);
    CHECK_NEAR(w.bodies[1].w.z - w.bodies[0].w.z, 2, 1e-3);  // b relative to a, about the axis
    CHECK_NEAR(w.bodies[1].w.z + w.bodies[0].w.z, 0, 1e-3);  // equal bodies: equal and opposite
    w.joints[static_cast<size_t>(motor)].motor_speed = -1;
    run(w, 60);
    CHECK_NEAR(w.bodies[1].w.z - w.bodies[0].w.z, -1, 1e-3);
}

TEST(a_motor_lifts_an_arm_only_if_its_torque_beats_the_weight) {
    auto angle_after = [](Real torque_over_needed) {
        World w;
        const int arm = w.add(Body::solid_box({0.4f, 0.05f, 0.05f}, 1, {0.4f, 0, 0}));
        Joint hinge = Joint::hinge(w.bodies, -1, arm, {0, 0, 0}, {0, 0, 1});
        hinge.enable_motor = true;
        hinge.motor_speed = 1;  // counter-clockwise seen from +z: lifts an arm pointing along +x
        hinge.max_motor = torque_over_needed * w.bodies[0].mass * kG * 0.4f;
        w.add_joint(hinge);
        run(w, 120);
        return w.joints[0].hinge_angle(w.bodies);
    };
    CHECK(angle_after(1.1f) > 0.02f);
    CHECK(angle_after(0.9f) < -0.02f);
}

TEST(hinge_limits_stop_a_falling_arm) {
    World w;
    const int arm = w.add(Body::solid_box({0.4f, 0.05f, 0.05f}, 1, {0.4f, 0, 0}));
    Joint hinge = Joint::hinge(w.bodies, -1, arm, {0, 0, 0}, {0, 0, 1});
    hinge.enable_limit = true;
    hinge.lower = -0.5f;
    hinge.upper = 0.25f;
    w.add_joint(hinge);
    Real lowest = 0;
    for (int i = 0; i < 240 * 3; ++i) {
        w.step(kDt);
        lowest = std::min(lowest, w.joints[0].hinge_angle(w.bodies));
    }
    CHECK(lowest > -0.5f - 0.01f);
    CHECK_NEAR(w.joints[0].hinge_angle(w.bodies), -0.5, 0.005);  // resting on the stop
    CHECK(w.bodies[0].w.length() < 0.01f);
    // The stop pushes one way only: the impulse holding the arm up is in the lower limit's slot.
    CHECK(w.joints[0].impulses[6] > 0);
    CHECK(w.joints[0].impulses[7] == 0);

    // Thrown upward, it stops at the upper limit instead.
    World up;
    up.gravity = {};
    up.add(Body::solid_box({0.4f, 0.05f, 0.05f}, 1, {0.4f, 0, 0}));
    up.add_joint(hinge);
    up.bodies[0].w = {0, 0, 2};
    up.bodies[0].vel = {0, 0.8f, 0};
    Real highest = 0;
    for (int i = 0; i < 240; ++i) {
        up.step(kDt);
        highest = std::max(highest, up.joints[0].hinge_angle(up.bodies));
    }
    CHECK(highest > 0.24f && highest < 0.25f + 0.01f);
}

// ---- in the world --------------------------------------------------------------------------------

TEST(jointed_bodies_do_not_collide_unless_asked) {
    auto contacts = [](bool collide_connected) {
        World w;
        w.gravity = {};
        const int a = w.add(Body::solid_box({0.3f, 0.3f, 0.3f}, 1, {0, 0, 0}));
        const int b = w.add(Body::solid_box({0.3f, 0.3f, 0.3f}, 1, {0.5f, 0, 0}));  // overlapping a
        Joint pin = Joint::ball(w.bodies, a, b, {0.25f, 0, 0});
        pin.collide_connected = collide_connected;
        w.add_joint(pin);
        w.step(kDt);
        return w.contacts().size();
    };
    CHECK(contacts(false) == 0);
    CHECK(contacts(true) == 1);
}

TEST(a_joint_makes_one_island_and_the_island_sleeps_and_wakes_together) {
    World w;
    const int upper = w.add(Body::solid_box({0.05f, 0.3f, 0.05f}, 1, {0, -0.3f, 0}));
    const int lower = w.add(Body::solid_box({0.05f, 0.3f, 0.05f}, 1, {0, -0.9f, 0}));
    w.add_joint(Joint::ball(w.bodies, -1, upper, {0, 0, 0}));
    Joint elbow = Joint::hinge(w.bodies, upper, lower, {0, -0.6f, 0}, {0, 0, 1});
    elbow.enable_motor = true;
    elbow.max_motor = 50;
    const int motor = w.add_joint(elbow);
    run(w, 240 * 2);
    CHECK(w.stats().islands == 1);
    CHECK(w.stats().awake_bodies == 0);  // hanging straight down at rest: asleep
    const Vec3 rest = w.bodies[1].pos;
    run(w, 240);
    CHECK(distance(w.bodies[1].pos, rest) == 0);  // and nothing creeps while it sleeps

    w.joints[static_cast<size_t>(motor)].motor_speed = 1;  // a command touches no velocity and no force
    run(w, 60);
    CHECK(w.stats().awake_bodies == 2);
    CHECK(std::fabs(w.joints[static_cast<size_t>(motor)].hinge_angle(w.bodies)) > 0.1f);
}

TEST(truncate_drops_joints_of_removed_bodies) {
    World w;
    const int a = w.add(Body::solid_sphere(0.1f, 1, {0, 0, 0}));
    const int b = w.add(Body::solid_sphere(0.1f, 1, {1, 0, 0}));
    w.add_joint(Joint::ball(w.bodies, -1, a, {0, 1, 0}));
    w.add_joint(Joint::distance(w.bodies, a, b, {0, 0, 0}, {1, 0, 0}));
    w.truncate(1);
    CHECK(w.joints.size() == 1);
    CHECK(w.joints[0].type == JointType::Ball);
    run(w, 10);  // and stepping afterwards is safe
}

TEST(a_ragdoll_arm_resting_on_the_floor_is_quiet) {
    // Joints and contacts in the same sweeps: a two-link arm pinned to the world, lying on the floor.
    World w;
    Body floor = Body::fixed_box({5, 0.5f, 5}, {0, -0.5f, 0});
    floor.restitution = 0;
    w.add(floor);
    const int upper = w.add(Body::solid_box({0.3f, 0.05f, 0.05f}, 1, {0.3f, 0.5f, 0}));
    const int lower = w.add(Body::solid_box({0.3f, 0.05f, 0.05f}, 1, {0.9f, 0.5f, 0}));
    const int shoulder = w.add_joint(Joint::ball(w.bodies, -1, upper, {0, 0.5f, 0}));
    const int elbow = w.add_joint(Joint::hinge(w.bodies, upper, lower, {0.6f, 0.5f, 0}, {0, 0, 1}));
    run(w, 240 * 4);
    CHECK(w.bodies[static_cast<size_t>(lower)].pos.y > 0.03f);  // on the floor, not through it
    CHECK(w.bodies[static_cast<size_t>(lower)].pos.y < 0.3f);   // and it did fall
    CHECK(gap(w, shoulder) < 2e-3f);
    CHECK(gap(w, elbow) < 2e-3f);
    CHECK(w.stats().awake_bodies == 0);  // came to rest and went to sleep
}
