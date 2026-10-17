#include "test.hpp"
#include <phys3d/world.hpp>

#include <phys/world.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

using namespace phys3d;

namespace {

constexpr Real kG = 9.81f;
constexpr Real kDt = 1.0f / 240.0f;

void run(World& w, int steps, Real dt = kDt) {
    for (int i = 0; i < steps; ++i) w.step(dt);
}
// Distance between a joint's two anchors: how far a pin joint has come apart.
Real gap_between(const std::vector<Body>& bodies, const Joint& j) { return distance(j.world_anchor_a(bodies), j.world_anchor_b(bodies)); }
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

TEST(pendulums_that_share_nothing_swing_bit_for_bit_as_they_would_alone) {
    // Three pendulums hung from the world, far apart: each must move exactly as it does in a world of its
    // own. The solver may take joints in any order (it colours them; see solving_order.hpp), and with no
    // moving body in common the order cannot change a single bit. It must also never mix up whose
    // velocities it is working on: all three are joints to the world, index -1, and so come one after
    // another with the same first body.
    auto hang = [](World& w, int kind, Real x) {
        const int bob = w.add(Body::solid_box({0.1f, 0.2f, 0.1f}, 1, {x + 0.6f, -0.8f, 0.3f}, Quat::from_axis_angle({1, 2, 3}, 0.4f)));
        w.bodies[static_cast<size_t>(bob)].w = {0.5f, -1, 2};
        const Vec3 pivot{x, 0, 0};
        if (kind == 0) w.add_joint(Joint::ball(w.bodies, -1, bob, pivot));
        if (kind == 1) w.add_joint(Joint::hinge(w.bodies, -1, bob, pivot, {0, 0, 1}));
        if (kind == 2) w.add_joint(Joint::distance(w.bodies, -1, bob, pivot, {x + 0.6f, -0.6f, 0.3f}));
        return bob;
    };
    World together;
    for (int kind = 0; kind < 3; ++kind) hang(together, kind, static_cast<Real>(kind) * 10);
    for (int i = 0; i < 300; ++i) together.step(kDt);
    for (int kind = 0; kind < 3; ++kind) {
        World alone;
        hang(alone, kind, static_cast<Real>(kind) * 10);
        for (int i = 0; i < 300; ++i) alone.step(kDt);
        const Body &a = alone.bodies[0], &b = together.bodies[static_cast<size_t>(kind)];
        CHECK(a.pos.x == b.pos.x && a.pos.y == b.pos.y && a.pos.z == b.pos.z);
        CHECK(a.q.x == b.q.x && a.q.y == b.q.y && a.q.z == b.q.z && a.q.w == b.q.w);
        CHECK(a.vel.x == b.vel.x && a.w.y == b.w.y && a.w.z == b.w.z);
        CHECK(a.vel.length() + a.w.length() > 0.1f);  // it is still moving
    }
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

TEST(one_correction_closes_an_anchor_gap_and_moves_no_centre_of_mass) {
    // The position pass corrects an anchor's three rows as one (step 21): it solves for the push that closes
    // the gap, the bodies' turning included, so one correction leaves only what the turning's curvature
    // adds, (gap / lever)^2 of the lever (measured: 0.3% of the gap; row by row it was 2.9%). And it pushes
    // the two bodies equally and oppositely, so their centre of mass stays where it was. (Not a hinge: its two
    // turning rows come after the anchor and open it again by turning the bodies back, row by row or not.
    // Solving them together with the anchor closes it, but holds colliding chains worse; see step 22.)
    for (const JointType type : {JointType::Ball, JointType::Weld}) {
        std::vector<Body> bodies;
        bodies.push_back(Body::solid_box({0.1f, 0.3f, 0.2f}, 1, {0, 0, 0}, Quat::from_axis_angle({1, 2, 3}, 0.4f)));
        bodies.push_back(Body::solid_box({0.2f, 0.1f, 0.3f}, 3, {0.5f, 0.2f, 0.1f}, Quat::from_axis_angle({3, 1, 2}, 0.7f)));
        const Vec3 anchor{0.25f, 0.1f, 0.05f};
        std::vector<Joint> joints{type == JointType::Ball ? Joint::ball(bodies, 0, 1, anchor) : Joint::weld(bodies, 0, 1, anchor)};
        bodies[1].pos += Vec3{0.003f, -0.004f, 0.002f};  // pulled apart, not turned: the anchor is all that is wrong
        const Real before = distance(joints[0].world_anchor_a(bodies), joints[0].world_anchor_b(bodies));
        const Vec3 centre = (bodies[0].pos * bodies[0].mass + bodies[1].pos * bodies[1].mass) / (bodies[0].mass + bodies[1].mass);
        SolverSettings once;
        once.joint_position_iterations = 1;
        solve_joint_positions(bodies, joints, once);
        CHECK(distance(joints[0].world_anchor_a(bodies), joints[0].world_anchor_b(bodies)) < before * 0.01f);
        const Vec3 now = (bodies[0].pos * bodies[0].mass + bodies[1].pos * bodies[1].mass) / (bodies[0].mass + bodies[1].mass);
        CHECK(distance(now, centre) < 1e-6f);
    }
}

TEST(a_badly_opened_anchor_is_eased_back_not_snapped) {
    // An anchor a metre open closes by 20 cm per correction, as each row of a joint does, whatever the joint.
    for (const JointType type : {JointType::Ball, JointType::Hinge, JointType::Weld}) {
        std::vector<Body> bodies;
        bodies.push_back(Body::solid_sphere(0.2f, 1, {0, 0, 0}));
        bodies.push_back(Body::solid_sphere(0.2f, 1, {0, 0, 0}));
        std::vector<Joint> joints{type == JointType::Ball    ? Joint::ball(bodies, 0, 1, {0, 0, 0})  // anchored at both centres
                                  : type == JointType::Hinge ? Joint::hinge(bodies, 0, 1, {0, 0, 0}, {0, 0, 1})
                                                             : Joint::weld(bodies, 0, 1, {0, 0, 0})};
        bodies[1].pos = {1, 0, 0};
        SolverSettings once;
        once.joint_position_iterations = 1;
        solve_joint_positions(bodies, joints, once);
        CHECK_NEAR(distance(bodies[0].pos, bodies[1].pos), 0.8, 1e-5);
    }
}

TEST(the_position_pass_brings_a_bent_cone_back_a_step_at_a_time) {
    // A ball joint whose cone allows 0.3 rad, bent to 0.6: one correction brings it back by at most the
    // 0.14 rad any turning row may take at once, and four of them, over the limit no more.
    std::vector<Body> bodies;
    bodies.push_back(Body::solid_box({0.1f, 0.4f, 0.1f}, 1, {0, 0, 0}));
    bodies.push_back(Body::solid_box({0.1f, 0.4f, 0.1f}, 1, {0, -0.8f, 0}));
    Joint shoulder = Joint::ball(bodies, 0, 1, {0, -0.4f, 0});
    shoulder.enable_cone = true;
    shoulder.cone_angle = 0.3f;
    const Vec3 pivot{0, -0.4f, 0};
    const Quat bend = Quat::from_axis_angle({0, 0, 1}, 0.6f);  // b swung 0.6 rad about the anchor
    bodies[1].q = bend;
    bodies[1].pos = pivot + rotate(bend, Vec3{0, -0.4f, 0});
    std::vector<Joint> joints{shoulder};
    CHECK_NEAR(joints[0].cone_tilt(bodies), 0.6, 1e-4);
    SolverSettings once;
    once.joint_position_iterations = 1;
    solve_joint_positions(bodies, joints, once);
    const Real tilt = joints[0].cone_tilt(bodies);
    CHECK(tilt < 0.6f - 0.1f && tilt > 0.6f - 0.15f);
    solve_joint_positions(bodies, joints, SolverSettings{});
    CHECK(joints[0].cone_tilt(bodies) < 0.3f + 1e-3f);
    CHECK(gap_between(bodies, joints[0]) < 1e-4f);
}

TEST(a_sleeping_joint_is_left_alone_even_when_open) {
    World w;
    const int bob = w.add(Body::solid_sphere(0.1f, 1, {0, -1, 0}));
    w.add_joint(Joint::ball(w.bodies, -1, bob, {0, 0, 0}));
    run(w, 240);
    CHECK(!w.bodies[0].awake);  // hanging still: asleep
    w.bodies[0].pos.x += 0.001f;  // moved without waking it
    const Vec3 moved = w.bodies[0].pos;
    run(w, 10);
    CHECK(!w.bodies[0].awake);
    CHECK(distance(w.bodies[0].pos, moved) == 0);
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

// ---- weld ----------------------------------------------------------------------------------------

namespace {

// The angle of the rotation that takes orientation `from` to `to`.
Real turn_between(const Quat& from, const Quat& to) { return rotation_angle(to * conjugate(from)); }

}  // namespace

TEST(weld_and_slider_rows_match_the_rate_their_errors_change) {
    const std::vector<Body> bodies = two_tumblers();
    // Six rows for a weld: three for the orientation, three for the point.
    check_rows_against_finite_differences(bodies, Joint::weld(bodies, 0, 1, {0.7f, 0.9f, 0.1f}), 6);
    check_rows_against_finite_differences(bodies, Joint::weld(bodies, -1, 1, {0.7f, 0.9f, 0.1f}), 6);

    // A slider: three for the orientation, two across the axis. Made, then b pushed well off the anchor,
    // both along the rail and off it: one term of the rows only exists when the anchors are apart.
    Joint slider = Joint::slider(bodies, 0, 1, {0.7f, 0.9f, 0.1f}, {1, -2, 0.5f});
    std::vector<Body> apart = bodies;
    apart[1].pos += Vec3{0.4f, -0.3f, 0.6f};
    check_rows_against_finite_differences(apart, slider, 5);
    // And with both ends of the rail passed, the two limit rows as well.
    slider.enable_limit = true;
    slider.lower = 5;
    slider.upper = -5;
    check_rows_against_finite_differences(apart, slider, 7);
    CHECK_NEAR(slider.slide(apart), dot(Vec3{0.4f, -0.3f, 0.6f}, Vec3{1, -2, 0.5f}.normalized()), 1e-5);
}

TEST(two_welded_bodies_tumble_as_one) {
    World w;
    w.gravity = {};
    w.allow_sleep = false;
    for (const Body& b : two_tumblers()) w.add(b);
    const Quat relative = conjugate(w.bodies[0].q) * w.bodies[1].q;
    const Vec3 offset = inv_rotate(w.bodies[0].q, w.bodies[1].pos - w.bodies[0].pos);  // b's centre, seen from a
    const int weld = w.add_joint(Joint::weld(w.bodies, 0, 1, {0.7f, 0.9f, 0.1f}));
    Vec3 p0, l0;
    for (const Body& b : w.bodies) {
        p0 += b.vel * b.mass;
        l0 += cross(b.pos, b.vel) * b.mass + b.angular_momentum();
    }
    Real worst_gap = 0, worst_turn = 0, worst_shift = 0;
    for (int i = 0; i < 240 * 4; ++i) {
        w.step(kDt);
        if (i < 5) continue;  // the first steps bring their different velocities into line
        worst_gap = std::max(worst_gap, gap(w, weld));
        worst_turn = std::max(worst_turn, turn_between(relative, conjugate(w.bodies[0].q) * w.bodies[1].q));
        worst_shift = std::max(worst_shift, distance(offset, inv_rotate(w.bodies[0].q, w.bodies[1].pos - w.bodies[0].pos)));
    }
    CHECK(worst_gap < 1e-3f);
    CHECK(worst_turn < 2e-3f);   // about a tenth of a degree
    CHECK(worst_shift < 2e-3f);  // b stays where it was in a's frame
    CHECK(distance(w.bodies[0].w, w.bodies[1].w) < 1e-3f);  // one body: one angular velocity
    Vec3 p1, l1;
    for (const Body& b : w.bodies) {
        p1 += b.vel * b.mass;
        l1 += cross(b.pos, b.vel) * b.mass + b.angular_momentum();
    }
    CHECK(distance(p0, p1) < 1e-3f * p0.length());
    CHECK(distance(l0, l1) < 0.03f * l0.length());
}

TEST(a_weld_does_not_care_which_of_its_two_names_an_orientation_goes_by) {
    // q and -q are the same orientation. Nothing may happen when a body's quaternion is swapped for its
    // negative: read naively, the "difference" from where it should be would be a full turn.
    World w;
    w.gravity = {};
    w.allow_sleep = false;
    for (Body b : two_tumblers()) {
        b.vel = {};
        b.w = {};
        w.add(b);
    }
    w.add_joint(Joint::weld(w.bodies, 0, 1, {0.7f, 0.9f, 0.1f}));
    w.add_joint(Joint::slider(w.bodies, -1, 0, {0, 0, 0}, {1, 0, 0}));
    const Vec3 p0 = w.bodies[0].pos, p1 = w.bodies[1].pos;
    const Quat q0 = w.bodies[0].q, q1 = w.bodies[1].q;
    w.bodies[1].q = q1 * Real(-1);
    w.bodies[0].q = q0 * Real(-1);
    JointSolver solver;
    solver.prepare(w.bodies, w.joints, kDt, SolverSettings{});
    for (const JointRow& r : solver.rows()) CHECK_NEAR(r.error, 0, 1e-5);
    run(w, 30);
    CHECK(distance(w.bodies[0].pos, p0) < 1e-5f && distance(w.bodies[1].pos, p1) < 1e-5f);
    CHECK(turn_between(q0, w.bodies[0].q) < 1e-4f && turn_between(q1, w.bodies[1].q) < 1e-4f);
    CHECK(w.bodies[1].w.length() < 1e-5f);
}

TEST(a_beam_welded_to_the_world_holds_its_weight_and_its_torque) {
    // A beam sticking out sideways from a weld at one end. The weld must supply an upward force equal to the
    // weight, and a torque equal to the weight times the distance to the beam's centre.
    World w;
    w.allow_sleep = false;
    const Real reach = 0.8f;
    const int beam = w.add(Body::solid_box({reach, 0.1f, 0.1f}, 2, {reach, 0, 0}));
    const int weld = w.add_joint(Joint::weld(w.bodies, -1, beam, {0, 0, 0}));
    run(w, 480);
    const Body& b = w.bodies[static_cast<size_t>(beam)];
    CHECK(distance(b.pos, Vec3{reach, 0, 0}) < 1e-3f);
    CHECK(turn_between(Quat{}, b.q) < 1e-3f);
    const Joint& j = w.joints[static_cast<size_t>(weld)];
    const Real weight = b.mass * kG * kDt;
    CHECK_NEAR(j.impulses[1], weight, 2e-3 * weight);                  // the point's y row: the weight
    CHECK_NEAR(j.impulses[5], weight * reach, 5e-3 * weight * reach);  // the orientation's z row: the torque
    CHECK_NEAR(j.impulses[0], 0, 1e-3 * weight);
    CHECK_NEAR(j.impulses[3], 0, 1e-3 * weight);
    CHECK_NEAR(j.impulses[4], 0, 1e-3 * weight);
}

// ---- slider --------------------------------------------------------------------------------------

TEST(a_body_on_a_slider_only_slides) {
    // A rail tilted 0.6 rad from straight down. Along it, gravity's share is g cos(0.6); across it, and in
    // turning, nothing is allowed, whatever the body is hit with.
    World w;
    w.allow_sleep = false;
    const Vec3 axis{std::sin(0.6f), -std::cos(0.6f), 0};
    const int block = w.add(Body::solid_box({0.2f, 0.3f, 0.25f}, 1, {1, 5, 2}, Quat::from_axis_angle({1, 2, 3}, 0.7f)));
    const Quat start = w.bodies[0].q;
    const int rail = w.add_joint(Joint::slider(w.bodies, -1, block, {1.1f, 5.2f, 1.9f}, axis));  // anchored off its centre
    w.bodies[0].vel = {2, 1, -3};  // mostly not along the rail
    w.bodies[0].w = {3, -2, 4};
    const Real along0 = dot(Vec3{2, 1, -3}, axis);
    Real worst_off_axis = 0, worst_turn = 0;
    const int steps = 240;
    for (int i = 0; i < steps; ++i) {
        w.step(kDt);
        const Joint& j = w.joints[static_cast<size_t>(rail)];
        const Vec3 d = j.world_anchor_b(w.bodies) - j.world_anchor_a(w.bodies);
        worst_off_axis = std::max(worst_off_axis, (d - axis * dot(d, axis)).length());
        worst_turn = std::max(worst_turn, turn_between(start, w.bodies[0].q));
    }
    CHECK(worst_off_axis < 2e-3f);
    CHECK(worst_turn < 3e-3f);
    CHECK(w.bodies[0].w.length() < 1e-3f);
    const Real t = static_cast<Real>(steps) * kDt, g_along = kG * std::cos(0.6f);
    CHECK_NEAR(dot(w.bodies[0].vel, axis), along0 + g_along * t, 0.005 * (along0 + g_along * t));
    CHECK_NEAR(w.joints[0].slide(w.bodies), along0 * t + 0.5 * g_along * t * t, 0.01 * (along0 * t + 0.5 * g_along * t * t));
    CHECK(cross(w.bodies[0].vel, axis).length() < 1e-3f);  // its velocity is along the rail and nothing else
}

TEST(a_sliders_limits_are_the_ends_of_its_rail) {
    World w;
    w.allow_sleep = false;
    const int block = w.add(Body::solid_box({0.2f, 0.2f, 0.2f}, 1, {0, 3, 0}));
    Joint rail = Joint::slider(w.bodies, -1, block, {0, 3, 0}, {0, -1, 0});  // positive is downward
    rail.enable_limit = true;
    rail.lower = -0.3f;
    rail.upper = 0.8f;
    w.add_joint(rail);
    Real deepest = 0;
    for (int i = 0; i < 480; ++i) {
        w.step(kDt);
        deepest = std::max(deepest, w.joints[0].slide(w.bodies));
    }
    CHECK(deepest < 0.8f + 0.005f);
    CHECK_NEAR(w.joints[0].slide(w.bodies), 0.8, 0.003);  // hanging on the bottom stop
    CHECK(w.bodies[0].vel.length() < 0.01f);
    CHECK(w.joints[0].impulses[7] > 0 && w.joints[0].impulses[6] == 0);  // held by the upper limit's row only
    CHECK_NEAR(w.joints[0].impulses[7], w.bodies[0].mass * kG * kDt, 1e-3 * w.bodies[0].mass * kG * kDt);

    // Thrown upward, it stops at the other end.
    w.bodies[0].vel = {0, 6, 0};
    Real highest = 1;
    for (int i = 0; i < 60; ++i) {
        w.step(kDt);
        highest = std::min(highest, w.joints[0].slide(w.bodies));
    }
    CHECK(highest < -0.29f && highest > -0.3f - 0.005f);
}

TEST(a_sliders_motor_pushes_with_a_limited_force) {
    // No gravity, from rest: at its force limit the motor gives v = F t / m.
    {
        World w;
        w.gravity = {};
        const int block = w.add(Body::solid_box({0.2f, 0.2f, 0.2f}, 1, {0, 0, 0}));
        Joint rail = Joint::slider(w.bodies, -1, block, {0.1f, 0.1f, 0}, {1, 0, 0});
        rail.enable_motor = true;
        rail.motor_speed = 1000;  // never reached
        rail.max_motor = 0.05f;
        w.add_joint(rail);
        run(w, 240);
        CHECK_NEAR(w.bodies[0].vel.x, 0.05 / w.bodies[0].mass, 0.002 * 0.05 / w.bodies[0].mass);
        // Strong enough, it simply holds the speed asked for, and reverses when told.
        w.joints[0].max_motor = 100;
        w.joints[0].motor_speed = 2;
        run(w, 30);
        CHECK_NEAR(w.bodies[0].vel.x, 2, 1e-3);
        w.joints[0].motor_speed = -1;
        run(w, 30);
        CHECK_NEAR(w.bodies[0].vel.x, -1, 1e-3);
    }
    // A lift: it raises its load only if its force beats the weight.
    auto climb = [](Real force_over_weight) {
        World w;
        const int cab = w.add(Body::solid_box({0.3f, 0.3f, 0.3f}, 1, {0, 2, 0}));
        Joint rail = Joint::slider(w.bodies, -1, cab, {0, 2, 0}, {0, 1, 0});
        rail.enable_motor = true;
        rail.motor_speed = 1;
        rail.max_motor = force_over_weight * w.bodies[0].mass * kG;
        w.add_joint(rail);
        run(w, 120);
        return w.joints[0].slide(w.bodies);
    };
    CHECK(climb(1.1f) > 0.02f);
    CHECK(climb(0.9f) < -0.02f);
}

TEST(a_sleeping_slider_wakes_when_its_motor_is_told_to_move) {
    World w;
    const int block = w.add(Body::solid_box({0.2f, 0.2f, 0.2f}, 1, {0, 3, 0}));
    Joint rail = Joint::slider(w.bodies, -1, block, {0, 3, 0}, {0, 1, 0});
    rail.enable_motor = true;  // a motor asked for no speed is a brake: it holds the block where it is
    rail.max_motor = 1000;
    w.add_joint(rail);
    run(w, 240);
    CHECK(w.stats().awake_bodies == 0);
    CHECK_NEAR(w.joints[0].slide(w.bodies), 0, 1e-3);
    w.joints[0].motor_speed = 0.5f;
    run(w, 240);  // one second
    CHECK(w.stats().awake_bodies == 1);
    CHECK_NEAR(w.joints[0].slide(w.bodies), 0.5, 0.01);
}

TEST(a_slider_matches_the_2d_engines_prismatic_joint) {
    // The same block on the same tilted rail with the same end stops, in both engines.
    const Real mass = 1.5f, tilt = 0.4f, inertia = 0.2f;
    World w3;
    w3.allow_sleep = false;
    w3.add(Body::dynamic(mass, Mat3::diagonal(inertia, inertia, inertia), {0, 4, 0}));
    Joint rail3 = Joint::slider(w3.bodies, -1, 0, {0, 4, 0}, {std::cos(tilt), -std::sin(tilt), 0});
    rail3.enable_limit = true;
    rail3.lower = -0.5f;
    rail3.upper = 1.2f;
    w3.add_joint(rail3);

    phys::World w2;
    w2.allow_sleep = false;
    phys::Body block(phys::Shape::make_polygon(phys::Polygon::box(0.2f, 0.2f)), {0, 4}, 0);
    block.mass = mass;
    block.inv_mass = 1 / mass;
    block.inertia = inertia;
    block.inv_inertia = 1 / inertia;
    w2.add(block);
    phys::Joint rail2 = phys::Joint::prismatic(w2.bodies, -1, 0, {0, 4}, {std::cos(tilt), -std::sin(tilt)});
    rail2.enable_limit = true;
    rail2.lower = -0.5f;
    rail2.upper = 1.2f;
    w2.add_joint(rail2);

    for (int i = 0; i < 300; ++i) {
        w3.step(kDt);
        w2.step(kDt);
        if (i == 100 || i == 299) {  // on the way down, and at rest on the stop
            CHECK_NEAR(w3.bodies[0].pos.x, w2.bodies[0].pos.x, 1e-3);
            CHECK_NEAR(w3.bodies[0].pos.y, w2.bodies[0].pos.y, 1e-3);
            CHECK_NEAR(w3.bodies[0].vel.x, w2.bodies[0].vel.x, 2e-3);
        }
    }
    CHECK_NEAR(w3.joints[0].slide(w3.bodies), 1.2, 0.003);
}

// ---- the cone on a ball joint ----------------------------------------------------------------------

TEST(cone_row_matches_the_rate_the_tilt_changes) {
    std::vector<Body> bodies = two_tumblers();
    Joint shoulder = Joint::ball(bodies, 0, 1, {0.7f, 0.9f, 0.1f});
    shoulder.limit_cone(bodies, {1, -2, 0.5f}, 0.3f);
    CHECK_NEAR(shoulder.cone_tilt(bodies), 0, 1e-4);
    // Bent 0.8 rad, well past the limit: the cone row is in force and reports how far past.
    bodies[1].q = Quat::from_axis_angle({0.3f, 1, -0.2f}, 0.8f) * bodies[1].q;
    CHECK(shoulder.cone_tilt(bodies) > 0.4f);
    check_rows_against_finite_differences(bodies, shoulder, 4);
    JointSolver solver;
    solver.prepare(bodies, std::span<Joint>(&shoulder, 1), kDt, SolverSettings{});
    CHECK_NEAR(solver.rows()[3].error, 0.3 - shoulder.cone_tilt(bodies), 1e-5);
}

TEST(a_cone_stops_a_ball_joint_bending_too_far_and_does_nothing_before_that) {
    auto swing = [](bool cone, Real push, Real* widest, Vec3* end) {
        World w;
        w.allow_sleep = false;
        const int arm = w.add(Body::solid_box({0.08f, 0.5f, 0.08f}, 1, {0, -0.5f, 0}));
        Joint shoulder = Joint::ball(w.bodies, -1, arm, {0, 0, 0});
        if (cone) shoulder.limit_cone(w.bodies, {0, -1, 0}, 0.5f);
        w.add_joint(shoulder);
        w.bodies[0].w = {0, 0, push};  // swung about the pivot, in the xy plane...
        w.bodies[0].vel = cross(Vec3{0, 0, push}, Vec3{0, -0.5f, 0});
        const Real e0 = total_energy(w), kinetic0 = w.bodies[0].kinetic_energy();
        *widest = 0;
        for (int i = 0; i < 240 * 3; ++i) {
            w.step(kDt);
            *widest = std::max(*widest, w.joints[0].cone_tilt(w.bodies));
            // (The integrator's own wobble is under 1% of the swing's energy; a stop that kicked would show.)
            CHECK(total_energy(w) <= e0 + 0.02f * kinetic0);
        }
        *end = w.bodies[0].pos;
        return gap(w, 0);
    };
    Real widest_free, widest_coned, widest_gentle_free, widest_gentle_coned;
    Vec3 end_free, end_coned, end_gentle_free, end_gentle_coned;
    // A hard push: without the cone it swings out to well over a radian; with it, to the cone and no further.
    swing(false, 5, &widest_free, &end_free);
    const Real coned_gap = swing(true, 5, &widest_coned, &end_coned);
    CHECK(widest_free > 1.0f);
    CHECK(widest_coned > 0.49f && widest_coned < 0.5f + 0.01f);
    CHECK(coned_gap < 2e-3f);
    // A gentle one never reaches the cone, and then the cone changes nothing at all.
    swing(false, 1, &widest_gentle_free, &end_gentle_free);
    swing(true, 1, &widest_gentle_coned, &end_gentle_coned);
    CHECK(widest_gentle_free < 0.45f && widest_gentle_free > 0.1f);
    CHECK(distance(end_gentle_free, end_gentle_coned) < 1e-5f);
}

TEST(a_cone_leaves_twisting_about_the_axis_free) {
    // Held against the cone by a sideways pull, and spinning about its own length: the cone must stop the
    // lean and not touch the spin.
    World w;
    w.gravity = {3, -9.81f, 0};
    w.allow_sleep = false;
    Body arm = Body::solid_capsule(0.5f, 0.1f, 1, {0, -0.6f, 0});
    arm.w = {0, 6, 0};
    w.add(arm);
    Joint shoulder = Joint::ball(w.bodies, -1, 0, {0, 0, 0});
    shoulder.limit_cone(w.bodies, {0, -1, 0}, 0.1f);  // the pull would lean it 0.3 rad; the cone allows 0.1
    w.add_joint(shoulder);
    run(w, 240 * 3);
    CHECK_NEAR(w.joints[0].cone_tilt(w.bodies), 0.1, 0.01);
    const Vec3 along = rotate(w.bodies[0].q, {0, 1, 0});
    CHECK_NEAR(std::fabs(dot(w.bodies[0].w, along)), 6, 0.1);
}

// ---- joint trees, solved exactly (step 23) ----------------------------------------------------------

namespace {

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Vec3 in_box(Real r) { return {range(-r, r), range(-r, r), range(-r, r)}; }
};

// The largest speed left in any rigid row that holds both ways (the rows the trees solve), after `sweeps`
// sweeps of the joint solver from the bodies' present velocities, with nothing carried over.
Real worst_two_way_speed(std::vector<Body>& bodies, std::vector<Joint>& joints, int sweeps) {
    JointSolver solver;
    SolverSettings settings;
    settings.warm_starting = false;
    solver.prepare(bodies, joints, kDt, settings);
    for (int i = 0; i < sweeps; ++i) solver.solve_velocity();
    Real worst = 0;
    for (const JointRow& r : solver.rows())
        if (r.lower == -std::numeric_limits<Real>::infinity() && r.upper == std::numeric_limits<Real>::infinity() && r.softness == 0)
            worst = std::max(worst, std::fabs(r.velocity() + r.bias));
    return worst;
}

// A chain of `links` boxes along x, joined end to end by joints of every rigid kind in turn, the first held to
// the world if `held`; each body given a velocity and spin of its own.
void make_chain(std::vector<Body>& bodies, std::vector<Joint>& joints, int links, bool held, Lcg& rng) {
    for (int i = 0; i < links; ++i) {
        Body b = Body::solid_box({0.3f, 0.08f, 0.1f}, 1 + rng.range(0, 2), {0.6f * static_cast<Real>(i) + 0.3f, 0, 0});
        b.vel = rng.in_box(2);
        b.w = rng.in_box(3);
        bodies.push_back(b);
    }
    for (int i = held ? 0 : 1; i < links; ++i) {
        const int a = i - 1;
        const Vec3 at{0.6f * static_cast<Real>(i), 0, 0};
        switch (i % 4) {
            case 0: joints.push_back(Joint::ball(bodies, a, i, at)); break;
            case 1: joints.push_back(Joint::hinge(bodies, a, i, at, Vec3{0, 1, 1}.normalized())); break;
            case 2: joints.push_back(Joint::weld(bodies, a, i, at)); break;
            default: joints.push_back(Joint::slider(bodies, a, i, at, {1, 0, 0})); break;
        }
    }
}

Vec3 momentum(const std::vector<Body>& bodies) {
    Vec3 p;
    for (const Body& b : bodies) p += b.vel * b.mass;
    return p;
}
Vec3 angular_momentum(const std::vector<Body>& bodies) {
    Vec3 l;
    for (const Body& b : bodies) l += cross(b.pos, b.vel * b.mass) + b.angular_momentum();
    return l;
}

}  // namespace

TEST(a_chain_of_every_rigid_joint_is_solved_exactly_in_one_sweep) {
    // Ball, hinge, weld and slider joints, twelve links, held to the world at one end or free: one sweep of
    // the joint solver leaves every two-way row still, where Gauss-Seidel, one row at a time, needs hundreds.
    for (const bool held : {true, false}) {
        Lcg rng;
        std::vector<Body> bodies;
        std::vector<Joint> joints;
        make_chain(bodies, joints, 12, held, rng);
        const Real exact = worst_two_way_speed(bodies, joints, 1);
        CHECK(exact < 1e-4f);
    }
}

TEST(a_branching_tree_of_joints_is_solved_exactly_in_one_sweep) {
    // A ragdoll's shape: a torso with a welded head and four limbs of two links each, at the shoulders and hips
    // ball joints, at the elbows and knees hinges. Held to nothing.
    Lcg rng;
    std::vector<Body> bodies;
    std::vector<Joint> joints;
    bodies.push_back(Body::solid_box({0.2f, 0.3f, 0.1f}, 1, {0, 0, 0}));
    bodies.push_back(Body::solid_sphere(0.12f, 1, {0, 0.45f, 0}));
    joints.push_back(Joint::weld(bodies, 0, 1, {0, 0.33f, 0}));
    const Vec3 roots[4] = {{0.25f, 0.25f, 0}, {-0.25f, 0.25f, 0}, {0.12f, -0.3f, 0}, {-0.12f, -0.3f, 0}};
    const Vec3 out[4] = {{1, 0, 0}, {-1, 0, 0}, {0, -1, 0}, {0, -1, 0}};
    for (int limb = 0; limb < 4; ++limb) {
        int previous = 0;
        for (int k = 0; k < 2; ++k) {
            const Vec3 start = roots[limb] + out[limb] * (0.3f * static_cast<Real>(k));
            const int link = static_cast<int>(bodies.size());
            bodies.push_back(Body::solid_capsule(0.1f, 0.05f, 1, start + out[limb] * 0.15f, Quat::from_axis_angle(cross({0, 1, 0}, out[limb]).length_sq() > 0 ? cross({0, 1, 0}, out[limb]).normalized() : Vec3{1, 0, 0}, kPi / 2)));
            joints.push_back(k == 0 ? Joint::ball(bodies, previous, link, start) : Joint::hinge(bodies, previous, link, start, {0, 0, 1}));
            previous = link;
        }
    }
    for (Body& b : bodies) b.vel = rng.in_box(2), b.w = rng.in_box(3);
    CHECK(worst_two_way_speed(bodies, joints, 1) < 1e-4f);
}

TEST(solving_a_free_tree_exactly_keeps_its_momentum) {
    // The impulses a joint applies are equal and opposite, so whatever they come to, they cannot change a free
    // structure's momentum, linear or angular.
    Lcg rng;
    std::vector<Body> bodies;
    std::vector<Joint> joints;
    make_chain(bodies, joints, 8, false, rng);
    const Vec3 p = momentum(bodies), l = angular_momentum(bodies);
    worst_two_way_speed(bodies, joints, 1);
    CHECK(distance(momentum(bodies), p) < 1e-4f * (p.length() + 1));
    CHECK(distance(angular_momentum(bodies), l) < 1e-4f * (l.length() + 1));
}

TEST(a_welded_cantilever_holds_its_shape) {
    // Ten links welded end to end and to the world, sticking out sideways and a little down. Solved one row at a
    // time, ten sweeps a step, it sagged (by 2 mm at the joints on average, 14 at worst); solved as a tree, it
    // does not move.
    World w;
    const Vec3 lean = Vec3{1, -0.3f, 0}.normalized();
    const Quat q = Quat::from_axis_angle(cross({0, -1, 0}, lean), std::acos(-lean.y));
    int previous = -1;
    for (int i = 0; i < 10; ++i) {
        const int link = w.add(Body::solid_box({0.08f, 0.35f, 0.08f}, 1, lean * (0.35f + 0.7f * static_cast<Real>(i)), q));
        w.add_joint(Joint::weld(w.bodies, previous, link, lean * (0.7f * static_cast<Real>(i))));
        previous = link;
    }
    std::vector<Vec3> start;
    for (const Body& b : w.bodies) start.push_back(b.pos);
    run(w, 240);
    for (size_t i = 0; i < w.bodies.size(); ++i) CHECK(distance(w.bodies[i].pos, start[i]) < 1e-4f);
}

TEST(loops_of_joints_are_solved_exactly_and_a_body_that_cannot_turn_is_left_to_sweeping) {
    // A ring of four bodies (a loop), a bridge of six held to the world at both ends (a loop through the world, and
    // straight, so that one of its cut rows is implied by the rest and is dropped), and a chain with a body that
    // cannot turn about one axis (whose mass matrix has no inverse). The loops are cut and put back exactly through
    // their trees (step 27): one sweep leaves them still. The third is left to Gauss-Seidel, which needs sweeping,
    // and must never turn that body about that axis.
    for (int kind = 0; kind < 3; ++kind) {
        Lcg rng;
        std::vector<Body> bodies;
        std::vector<Joint> joints;
        const int n = kind == 0 ? 4 : 6;
        for (int i = 0; i < n; ++i) {
            const Real angle = 2 * kPi * static_cast<Real>(i) / static_cast<Real>(n);
            const Vec3 at = kind == 0 ? Vec3{std::cos(angle), std::sin(angle), 0} : Vec3{0.6f * static_cast<Real>(i) + 0.3f, 0, 0};
            bodies.push_back(Body::solid_box({0.2f, 0.08f, 0.1f}, 1, at));
            bodies.back().vel = rng.in_box(1);
            bodies.back().w = rng.in_box(1);
        }
        if (kind == 2) bodies[2].inv_inertia_body = Mat3::diagonal(1, 1, 0);  // (cannot turn about z)
        for (int i = 0; i < n; ++i) {
            if (kind == 0) joints.push_back(Joint::ball(bodies, i, (i + 1) % n, (bodies[static_cast<size_t>(i)].pos + bodies[static_cast<size_t>((i + 1) % n)].pos) * 0.5f));
            else joints.push_back(Joint::ball(bodies, i - 1, i, {0.6f * static_cast<Real>(i), 0, 0}));
        }
        if (kind == 1) joints.push_back(Joint::ball(bodies, n - 1, -1, {0.6f * static_cast<Real>(n), 0, 0}));
        const Real locked_spin = bodies[2].w.z;
        std::vector<Body> again = bodies;
        const Real after_one = worst_two_way_speed(bodies, joints, 1);
        const Real after_many = worst_two_way_speed(again, joints, 400);
        CHECK(std::isfinite(after_one) && std::isfinite(after_many));
        CHECK(after_many < 1e-3f);
        if (kind < 2) CHECK(after_one < 1e-4f);  // exactly, at once
        else CHECK(after_one > 1e-2f);           // (by sweeping)
        if (kind == 2) CHECK(again[2].w.z == locked_spin);  // no impulse turned it about its locked axis
    }
}

// ---- joint trees in the position pass (step 26) ------------------------------------------------------

namespace {

Real worst_gap(const std::vector<Body>& bodies, const std::vector<Joint>& joints) {
    Real worst = 0;
    for (const Joint& j : joints) worst = std::max(worst, gap_between(bodies, j));
    return worst;
}

// The largest position error left in any rigid row that holds both ways (a slider's anchors may part along its
// rail; its rows say what it holds).
Real worst_error(std::vector<Body>& bodies, std::vector<Joint>& joints) {
    JointSolver solver;
    solver.prepare(bodies, joints, kDt, SolverSettings{}, /*for_velocity=*/false);
    Real worst = 0;
    for (const JointRow& r : solver.rows())
        if (r.lower == -std::numeric_limits<Real>::infinity() && r.upper == std::numeric_limits<Real>::infinity()) worst = std::max(worst, std::fabs(r.error));
    return worst;
}

}  // namespace

TEST(one_position_pass_closes_a_whole_tree) {
    // A chain of every rigid joint, its links knocked a few millimetres and a few milliradians out of place: the
    // position pass corrects all of a tree at once, from one measurement, so one pass leaves only what the
    // turning's curvature adds, micrometres. (One joint at a time it left a tenth of the gaps: each joint's
    // correction opened its neighbours.)
    Lcg rng;
    std::vector<Body> bodies;
    std::vector<Joint> joints;
    make_chain(bodies, joints, 12, true, rng);
    for (Body& b : bodies) {
        b.pos += rng.in_box(0.004f);
        b.q = (Quat::from_axis_angle(Vec3{rng.range(-1, 1), rng.range(-1, 1), 1}.normalized(), rng.range(-0.004f, 0.004f)) * b.q).normalized();
    }
    const Real before = worst_error(bodies, joints);
    SolverSettings once;
    once.joint_position_iterations = 1;
    solve_joint_positions(bodies, joints, once);
    CHECK(before > 5e-3f);
    CHECK(worst_error(bodies, joints) < 5e-5f);  // (measured: 1.2 * 10^-5)
}

TEST(a_tree_knocked_far_out_of_place_gets_more_passes) {
    // Knocked by centimetres, one correction leaves what the curvature adds, a millimetre or so; a tree whose
    // errors were that large is corrected again in the next pass (and only such a tree), until they are not.
    Lcg rng;
    std::vector<Body> bodies;
    std::vector<Joint> joints;
    make_chain(bodies, joints, 8, true, rng);
    for (Body& b : bodies) b.pos += rng.in_box(0.05f);
    SolverSettings once, four;
    once.joint_position_iterations = 1;
    four.joint_position_iterations = 4;
    std::vector<Body> again = bodies;
    solve_joint_positions(bodies, joints, once);
    solve_joint_positions(again, joints, four);
    CHECK(worst_error(bodies, joints) > 1e-4f);  // one pass is not enough for this
    CHECK(worst_error(again, joints) < 1e-5f);
}

TEST(a_loop_of_joints_is_closed_in_the_position_pass_like_a_tree) {
    // A ring of four bodies joined in a loop, one of them pulled 2.4 cm out: the position pass cuts the loop and puts
    // the cut back through the tree, so one pass leaves what the turning's curvature adds, and a second next to
    // nothing. (Row by row it took twenty passes to get to 4 mm.)
    std::vector<Body> bodies;
    std::vector<Joint> joints;
    for (int i = 0; i < 4; ++i) {
        const Real angle = kPi / 2 * static_cast<Real>(i);
        bodies.push_back(Body::solid_box({0.2f, 0.08f, 0.1f}, 1, {std::cos(angle), std::sin(angle), 0}));
    }
    for (int i = 0; i < 4; ++i)
        joints.push_back(Joint::ball(bodies, i, (i + 1) % 4, (bodies[static_cast<size_t>(i)].pos + bodies[static_cast<size_t>((i + 1) % 4)].pos) * 0.5f));
    bodies[1].pos += Vec3{0.01f, 0.02f, -0.01f};
    CHECK(worst_gap(bodies, joints) > 0.02f);
    SolverSettings once;
    once.joint_position_iterations = 1;
    solve_joint_positions(bodies, joints, once);
    CHECK(worst_gap(bodies, joints) < 2e-4f);  // (measured: 5 * 10^-5)
    solve_joint_positions(bodies, joints, once);
    CHECK(worst_gap(bodies, joints) < 1e-6f);
}

TEST(a_rope_bridge_held_at_both_ends_hangs_together) {
    // Twelve planks joined by ball joints, the ends joined to the world: a loop through the world. Started as a V
    // with slack (six planks 30 degrees down, six up), it swings and settles into a sag; its joints stay closed.
    World w;
    const Real tilt = kPi / 6, length = 0.45f;
    int previous = -1;
    Vec3 at{0, 5, 0};
    for (int i = 0; i < 12; ++i) {
        const Real angle = i < 6 ? -tilt : tilt;
        const Vec3 along{std::cos(angle), std::sin(angle), 0};
        const int plank = w.add(Body::solid_box({length / 2, 0.03f, 0.4f}, 1, at + along * (length / 2), Quat::from_axis_angle({0, 0, 1}, angle)));
        w.add_joint(Joint::ball(w.bodies, previous, plank, at));
        previous = plank;
        at += along * length;
    }
    w.add_joint(Joint::ball(w.bodies, previous, -1, at));
    Real worst = 0;
    for (int step = 0; step < 600; ++step) {
        w.step(kDt);
        worst = std::max(worst, worst_gap(w.bodies, w.joints));
    }
    CHECK(w.bodies[5].pos.y < 4.0f);  // it sags
    CHECK(worst < 1e-4f);             // and holds together
}

TEST(a_resting_rope_bridge_stores_exactly_the_impulse_that_carries_its_weight) {
    // The impulses kept for next step's warm start, a loop's included (one end's joint is a tree's, the other the
    // cut, put back through it): settled, the two ends between them carry the bridge's weight each step.
    World w;
    w.allow_sleep = false;
    const Real tilt = kPi / 6, length = 0.45f;
    int previous = -1;
    Vec3 at{0, 5, 0};
    for (int i = 0; i < 10; ++i) {
        const Real angle = i < 5 ? -tilt : tilt;
        const Vec3 along{std::cos(angle), std::sin(angle), 0};
        const int plank = w.add(Body::solid_box({length / 2, 0.03f, 0.4f}, 1, at + along * (length / 2), Quat::from_axis_angle({0, 0, 1}, angle)));
        w.add_joint(Joint::ball(w.bodies, previous, plank, at));
        previous = plank;
        at += along * length;
    }
    const int last = w.add_joint(Joint::ball(w.bodies, previous, -1, at));
    for (int step = 0; step < 960; ++step) {  // (settled by draining its motion by hand)
        w.step(kDt);
        for (Body& b : w.bodies) b.vel *= 0.97f, b.w *= 0.97f;
    }
    run(w, 1);
    Real weight = 0;
    for (const Body& b : w.bodies) weight += b.mass * kG * kDt;
    // The first joint pushes its plank (b) up by +y; the last pushes its plank (a) by -y.
    const Real held = w.joints[0].impulses[1] - w.joints[static_cast<size_t>(last)].impulses[1];
    CHECK_NEAR(held, weight, 2e-3 * weight);
}

TEST(the_impulses_kept_are_the_impulses_that_were_applied) {
    // Each row's accumulated impulse is what the next step starts from, so it must be what the row really did. From
    // a cold start (nothing carried over), one sweep: Gauss-Seidel, then the trees and the loops exactly. Then each
    // body's velocity must be where it began plus what every row's kept impulse gives it, M^-1 J^T lambda, for the
    // rows of trees, of loops' cuts, and of everything else alike.
    for (int kind = 0; kind < 3; ++kind) {
        Lcg rng;
        std::vector<Body> bodies;
        std::vector<Joint> joints;
        if (kind == 0) {
            make_chain(bodies, joints, 8, true, rng);  // a tree
        } else {
            const int n = kind == 1 ? 4 : 6;  // a ring; a bridge held at both ends
            for (int i = 0; i < n; ++i) {
                const Real angle = 2 * kPi * static_cast<Real>(i) / static_cast<Real>(n);
                bodies.push_back(Body::solid_box({0.2f, 0.08f, 0.1f}, 1, kind == 1 ? Vec3{std::cos(angle), std::sin(angle), 0} : Vec3{0.6f * static_cast<Real>(i) + 0.3f, 0.1f * static_cast<Real>(i % 2), 0}));
                bodies.back().vel = rng.in_box(1);
                bodies.back().w = rng.in_box(1);
            }
            for (int i = 0; i < n; ++i) {
                if (kind == 1) joints.push_back(Joint::ball(bodies, i, (i + 1) % n, (bodies[static_cast<size_t>(i)].pos + bodies[static_cast<size_t>((i + 1) % n)].pos) * 0.5f));
                else joints.push_back(Joint::ball(bodies, i - 1, i, i == 0 ? Vec3{0, 0, 0} : (bodies[static_cast<size_t>(i - 1)].pos + bodies[static_cast<size_t>(i)].pos) * 0.5f));
            }
            if (kind == 2) joints.push_back(Joint::ball(bodies, n - 1, -1, bodies.back().pos + Vec3{0.3f, 0, 0}));
        }
        std::vector<Vec3> v0, w0;
        for (const Body& b : bodies) v0.push_back(b.vel), w0.push_back(b.w);
        JointSolver solver;
        SolverSettings cold;
        cold.warm_starting = false;
        solver.prepare(bodies, joints, kDt, cold);
        solver.solve_velocity(/*exactly=*/true);
        std::vector<Vec3> v = v0, w = w0;
        for (const JointRow& r : solver.rows()) {
            if (r.a->type == BodyType::Dynamic) {
                const size_t i = static_cast<size_t>(r.a - bodies.data());
                v[i] += r.lin_a * (r.a->inv_mass * r.impulse), w[i] += r.inv_i_ang_a * r.impulse;
            }
            if (r.b->type == BodyType::Dynamic) {
                const size_t i = static_cast<size_t>(r.b - bodies.data());
                v[i] += r.lin_b * (r.b->inv_mass * r.impulse), w[i] += r.inv_i_ang_b * r.impulse;
            }
        }
        Real worst = 0;
        for (size_t i = 0; i < bodies.size(); ++i) worst = std::max({worst, distance(v[i], bodies[i].vel), distance(w[i], bodies[i].w)});
        CHECK(worst < 1e-4f);
    }
}

TEST(a_joint_given_twice_is_no_trouble) {
    // Two identical ball joints between the same two bodies: the second closes a loop, and every one of its rows is
    // implied by the first's, so all are dropped (given no impulse). The pair is still solved exactly.
    Lcg rng;
    std::vector<Body> bodies;
    bodies.push_back(Body::solid_box({0.2f, 0.1f, 0.1f}, 1, {0, 0, 0}));
    bodies.push_back(Body::solid_box({0.2f, 0.1f, 0.1f}, 1, {0.4f, 0, 0}));
    for (Body& b : bodies) b.vel = rng.in_box(1), b.w = rng.in_box(2);
    std::vector<Joint> joints{Joint::ball(bodies, 0, 1, {0.2f, 0, 0}), Joint::ball(bodies, 0, 1, {0.2f, 0, 0})};
    const Real left = worst_two_way_speed(bodies, joints, 1);
    CHECK(std::isfinite(left) && left < 1e-4f);
    CHECK(std::isfinite(bodies[0].vel.x) && bodies[0].vel.length() < 10);
}
