#include "test.hpp"
#include <phys/world.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

using namespace phys;

namespace {

constexpr double kG = 9.81;
constexpr double kPiD = 3.14159265358979323846;

Body disc(Real radius, Vec2 pos, Real angle = 0) { return Body(Shape::make_circle(radius), pos, angle); }
Body plank(Real hx, Real hy, Vec2 pos, Real angle = 0) {
    return Body(Shape::make_polygon(Polygon::box(hx, hy)), pos, angle);
}

void run(World& w, int steps, Real dt) {
    for (int i = 0; i < steps; ++i) w.step(dt);
}

// Separation of a joint's two anchors: how badly the joint is currently violated.
double anchor_error(const World& w, const Joint& j) {
    return double((j.world_anchor_b(w.bodies) - j.world_anchor_a(w.bodies)).length());
}

}  // namespace

// ---- the Jacobians themselves --------------------------------------------------------------------

// Each row claims "J v is the rate of change of the constrained quantity". Check that against a
// numerical derivative: move the bodies a hair along their velocities and see how the constraint
// values really change. A wrong sign or lever arm in any row shows up here immediately.
namespace {

using ConstraintFn = std::function<std::vector<double>(const std::vector<Body>&)>;

void check_rows_match_finite_differences(const Joint& joint, std::vector<Body> bodies, const ConstraintFn& values,
                                         size_t expected_rows) {
    std::vector<Joint> joints = {joint};
    JointSolver solver;
    solver.prepare(bodies, joints, 1.0f / 60.0f, SolverSettings{});
    CHECK(solver.rows().size() == expected_rows);

    auto moved = [&](Real h) {
        std::vector<Body> b = bodies;
        for (Body& body : b) {
            body.pos += body.vel * h;
            body.set_angle(body.angle + body.w * h);
        }
        return values(b);
    };
    const Real h = 1e-3f;
    const std::vector<double> plus = moved(h), minus = moved(-h);
    CHECK(plus.size() == solver.rows().size());
    for (size_t i = 0; i < solver.rows().size() && i < plus.size(); ++i) {
        const JointRow& r = solver.rows()[i];
        const double jv = double(dot(r.lin_a, r.a->vel)) + double(r.ang_a * r.a->w) + double(dot(r.lin_b, r.b->vel)) +
                          double(r.ang_b * r.b->w);
        const double numeric = (plus[i] - minus[i]) / (2.0 * double(h));
        CHECK_NEAR(jv, numeric, 8e-3 * (1.0 + std::fabs(numeric)));
    }
}

std::vector<Body> two_moving_bodies() {
    std::vector<Body> b = {plank(0.5f, 0.3f, {1.0f, 0.5f}, 0.5f), plank(0.4f, 0.4f, {3.2f, 1.3f}, -0.8f)};
    b[0].vel = {0.7f, -0.4f};
    b[0].w = 0.9f;
    b[1].vel = {-0.3f, 0.8f};
    b[1].w = -1.3f;
    return b;
}

}  // namespace

TEST(revolute_rows_are_the_derivative_of_the_anchor_separation) {
    auto bodies = two_moving_bodies();
    Joint j = Joint::revolute(bodies, 0, 1, {2.0f, 0.9f});
    j.enable_motor = true;
    j.enable_limit = true;
    j.lower = -2;
    j.upper = 2;
    check_rows_match_finite_differences(
        j, bodies,
        [&](const std::vector<Body>& b) {
            const Vec2 d = j.world_anchor_b(b) - j.world_anchor_a(b);
            const double rel = double(b[1].angle) - double(b[0].angle);
            return std::vector<double>{double(d.x), double(d.y), rel, rel, -rel};  // x, y, motor, lower, upper
        },
        5);
}

TEST(distance_rows_are_the_derivative_of_the_anchor_distance) {
    auto bodies = two_moving_bodies();
    Joint between = Joint::distance(bodies, 0, 1, {1.3f, 0.8f}, {2.9f, 1.0f});
    check_rows_match_finite_differences(
        between, bodies,
        [&](const std::vector<Body>& b) {
            return std::vector<double>{double((between.world_anchor_b(b) - between.world_anchor_a(b)).length())};
        },
        1);

    Joint to_world = Joint::distance(bodies, -1, 1, {0.0f, 4.0f}, {3.2f, 1.6f});  // a = the fixed world
    check_rows_match_finite_differences(
        to_world, bodies,
        [&](const std::vector<Body>& b) {
            return std::vector<double>{double((to_world.world_anchor_b(b) - to_world.world_anchor_a(b)).length())};
        },
        1);
}

TEST(prismatic_rows_are_the_derivative_of_their_constraints) {
    auto bodies = two_moving_bodies();
    Joint j = Joint::prismatic(bodies, 0, 1, {2.0f, 0.9f}, {1.0f, 0.4f});
    // Slide body b along the rail and off it, so the anchors no longer coincide. The rail's own turning
    // (it is fixed in body a) contributes a Jacobian term proportional to that offset, which would be
    // invisible with zero offset.
    bodies[1].pos += Vec2{0.5f, 0.35f};
    j.enable_motor = true;
    j.enable_limit = true;
    j.lower = -3;
    j.upper = 3;
    check_rows_match_finite_differences(
        j, bodies,
        [&](const std::vector<Body>& b) {
            const Vec2 d = j.world_anchor_b(b) - j.world_anchor_a(b);
            const Vec2 axis = rotate(b[0].q, j.local_axis_a);
            const double t = double(dot(axis, d));
            const double off_rail = double(dot(perp(axis), d));
            const double twist = double(b[1].angle) - double(b[0].angle);
            return std::vector<double>{off_rail, twist, t, t, -t};  // perp, no-rotation, motor, lower, upper
        },
        5);
}

TEST(mouse_rows_are_the_derivative_of_the_grab_point_error) {
    auto bodies = two_moving_bodies();
    Joint j = Joint::mouse(bodies, 1, {3.4f, 1.5f});
    j.target = {5.0f, 4.0f};
    check_rows_match_finite_differences(
        j, bodies,
        [&](const std::vector<Body>& b) {
            const Vec2 e = j.world_anchor_a(b) - j.target;
            return std::vector<double>{double(e.x), double(e.y)};
        },
        2);
}

// ---- revolute -------------------------------------------------------------------------------------

// A disc of radius r hanging by a pin L above its centre is a physical pendulum with
// I_pivot = m r^2/2 + m L^2. For amplitude theta0 the period is T0 (1 + theta0^2/16 + ...).
TEST(pinned_disc_swings_with_the_physical_pendulum_period) {
    const Real r = 0.2f, L = 1.5f, theta0 = 0.3f, dt = 1.0f / 240.0f;
    World w;
    const Vec2 pivot{0, 5};
    w.add(disc(r, {pivot.x - L * std::sin(theta0), pivot.y - L * std::cos(theta0)}));
    w.add_joint(Joint::revolute(w.bodies, -1, 0, pivot));

    const double t0 = 2 * kPiD * std::sqrt((0.5 * r * r + double(L) * L) / (kG * L));
    const double expected = t0 * (1 + 0.3 * 0.3 / 16.0 + 11.0 * std::pow(0.3, 4) / 3072.0);

    // Upward zero crossings of x - pivot.x are one period apart.
    std::vector<double> crossings;
    double max_error = 0;
    Real prev = w.bodies[0].pos.x - pivot.x;
    for (int i = 1; i <= 240 * 12; ++i) {
        w.step(dt);
        const Real x = w.bodies[0].pos.x - pivot.x;
        if (prev < 0 && x >= 0) crossings.push_back((double(i) - double(x) / double(x - prev)) * double(dt));
        prev = x;
        max_error = std::max(max_error, anchor_error(w, w.joints[0]));
    }
    CHECK(crossings.size() >= 4);
    if (crossings.size() >= 4) {
        const double period = (crossings.back() - crossings.front()) / double(crossings.size() - 1);
        CHECK_NEAR(period, expected, 0.01 * expected);
    }
    CHECK(max_error < 0.005);  // the hinge never comes apart by more than 5 mm
}

TEST(a_hanging_chain_stays_connected_while_swinging) {
    auto build = [](World& w, int links) {
        const Vec2 top{0, 8};
        for (int i = 0; i < links; ++i) w.add(plank(0.05f, 0.25f, {top.x, top.y - 0.25f - 0.5f * static_cast<Real>(i)}));
        w.add_joint(Joint::revolute(w.bodies, -1, 0, top));
        for (int i = 1; i < links; ++i)
            w.add_joint(Joint::revolute(w.bodies, i - 1, i, {top.x, top.y - 0.5f * static_cast<Real>(i)}));
        w.bodies.back().vel = {4.0f, 0};  // knock the free end sideways
    };
    auto worst_error = [&](World& w) {
        double worst = 0;
        for (int i = 0; i < 600; ++i) {
            w.step(1.0f / 120.0f);
            for (const Joint& j : w.joints) worst = std::max(worst, anchor_error(w, j));
        }
        return worst;
    };
    World w;
    build(w, 8);
    CHECK(worst_error(w) < 0.03);  // 8 hinges, a hard knock: no joint more than 3 cm apart

    // Warm starting carries each joint's impulse over, so with few sweeps it holds the chain together
    // markedly better.
    World warm, cold;
    warm.solver.iterations = cold.solver.iterations = 3;
    cold.solver.warm_starting = false;
    build(warm, 8);
    build(cold, 8);
    CHECK(worst_error(warm) < worst_error(cold));
}

TEST(revolute_limit_stops_the_swing) {
    World w;
    const Vec2 pivot{0, 5};
    w.add(plank(0.1f, 0.5f, {pivot.x, pivot.y - 0.5f}));  // a rod hanging from the pivot
    Joint j = Joint::revolute(w.bodies, -1, 0, pivot);
    j.enable_limit = true;
    j.lower = -0.5f;
    j.upper = 0.5f;
    w.add_joint(j);
    // Swung hard counter-clockwise about the pivot: w = 6 rad/s, with the centre of mass (0.5 m below
    // the pivot) moving at w x r = 3 m/s to match, as the pin demands.
    w.bodies[0].w = 6;
    w.bodies[0].vel = {3, 0};

    double max_angle = -10, min_angle = 10;
    for (int i = 0; i < 480; ++i) {
        w.step(1.0f / 120.0f);
        max_angle = std::max(max_angle, double(w.bodies[0].angle));
        min_angle = std::min(min_angle, double(w.bodies[0].angle));
    }
    CHECK(max_angle > 0.45);   // it really reached the limit...
    CHECK(max_angle < 0.56);   // ...and was stopped there (a little give is Baumgarte's)
    CHECK(min_angle > -0.56);  // and the swing back did not break the other limit
}

// An arm pinned at one end falls onto its lower limit and must come to rest THERE, pinned. The position
// pass used to correct all of a joint's rows from one measurement: closing the pin turned the arm off the
// limit, turning it back opened the pin, and the two overshot against each other for ever. The arm sat a
// few degrees above its stop with the pin pulled apart, and since it was "not at the limit" nothing ever
// stopped its speed building up. Over half of these cases failed.
TEST(an_arm_dropped_onto_a_revolute_limit_rests_exactly_on_it) {
    for (Real hz : {60.0f, 120.0f, 240.0f})
        for (Real half_length : {0.2f, 0.8f, 1.5f})
            for (Real half_thickness : {0.02f, 0.2f})
                for (Real lower : {-0.2f, -0.5f, -1.2f}) {
                    World w;
                    w.allow_sleep = false;
                    w.add(Body(Shape::make_polygon(Polygon::box(half_length, half_thickness)), {half_length, 0}, 0));
                    Joint hinge = Joint::revolute(w.bodies, -1, 0, {0, 0});
                    hinge.enable_limit = true;
                    hinge.lower = lower;
                    hinge.upper = 0.25f;
                    w.add_joint(hinge);
                    run(w, static_cast<int>(hz) * 4, 1 / hz);
                    CHECK_NEAR(w.bodies[0].angle, lower, 1e-3);
                    CHECK(std::fabs(w.bodies[0].w) < 1e-3f);
                    CHECK(distance(w.joints[0].world_anchor_a(w.bodies), w.joints[0].world_anchor_b(w.bodies)) < 1e-4f);
                }
}

TEST(revolute_motor_spins_to_speed_and_respects_its_torque_limit) {
    // A free wheel pinned at its centre. With plenty of torque it reaches the target speed...
    {
        World w;
        w.gravity = {};
        w.add(disc(0.5f, {0, 0}));
        Joint j = Joint::revolute(w.bodies, -1, 0, {0, 0});
        j.enable_motor = true;
        j.motor_speed = 3;
        j.max_motor = 1000;
        w.add_joint(j);
        run(w, 120, 1.0f / 120.0f);
        CHECK_NEAR(w.bodies[0].w, 3, 0.01);
    }
    // ...with little torque it accelerates at tau / I instead. (Disc: I = m r^2 / 2.)
    {
        World w;
        w.gravity = {};
        w.add(disc(0.5f, {0, 0}));
        const double inertia = double(w.bodies[0].inertia);
        Joint j = Joint::revolute(w.bodies, -1, 0, {0, 0});
        j.enable_motor = true;
        j.motor_speed = 30;
        j.max_motor = 0.5f;
        w.add_joint(j);
        run(w, 48, 1.0f / 120.0f);  // 0.4 s
        CHECK_NEAR(w.bodies[0].w, 0.5 / inertia * 0.4, 0.03 * 0.5 / inertia * 0.4);
    }
}

// ---- distance -------------------------------------------------------------------------------------

// A dumbbell: two bodies joined by a rigid rod, thrown so it spins. The rod must hold its length, and
// since the joint only exerts internal forces the total momentum cannot change.
TEST(rigid_distance_joint_holds_length_and_conserves_momentum) {
    World w;
    w.gravity = {};
    w.add(disc(0.2f, {0, 0}));
    w.add(disc(0.2f, {2, 0}));
    w.add_joint(Joint::distance(w.bodies, 0, 1, {0, 0}, {2, 0}));
    w.bodies[0].vel = {0.5f, -1};
    w.bodies[1].vel = {0.5f, 1};
    const double momentum_x = double(w.bodies[0].mass) * 0.5 + double(w.bodies[1].mass) * 0.5;

    double worst = 0;
    for (int i = 0; i < 600; ++i) {
        w.step(1.0f / 120.0f);
        worst = std::max(worst, std::fabs(double(distance(w.bodies[0].pos, w.bodies[1].pos)) - 2.0));
    }
    CHECK(worst < 0.01);
    const Body &a = w.bodies[0], &b = w.bodies[1];
    CHECK_NEAR(double(a.mass) * double(a.vel.x) + double(b.mass) * double(b.vel.x), momentum_x, 1e-4);
    CHECK_NEAR(double(a.mass) * double(a.vel.y) + double(b.mass) * double(b.vel.y), 0, 1e-4);
    // It is rotating about its centre of mass: both ends keep moving.
    CHECK(a.vel.length() > 0.5);
}

namespace {

// A disc on a spring: anchored to the world at the origin, stretched 0.3 m beyond its 2 m rest length.
// Returns the x positions sampled every step.
std::vector<double> spring_run(Real hz, Real zeta, int steps) {
    World w;
    w.gravity = {};
    w.add(disc(0.2f, {2.3f, 0}));
    Joint j = Joint::distance(w.bodies, -1, 0, {0, 0}, {2.3f, 0});
    j.length = 2;
    j.frequency_hz = hz;
    j.damping_ratio = zeta;
    w.add_joint(j);
    std::vector<double> xs;
    for (int i = 0; i < steps; ++i) {
        w.step(1.0f / 240.0f);
        xs.push_back(double(w.bodies[0].pos.x));
    }
    return xs;
}

}  // namespace

TEST(soft_distance_joint_oscillates_at_its_frequency) {
    const auto xs = spring_run(2.0f, 0.0f, 240 * 3);
    // Downward zero crossings of (x - rest length): consecutive ones are one period (0.5 s) apart.
    std::vector<double> crossings;
    for (size_t i = 1; i < xs.size(); ++i)
        if (xs[i - 1] > 2 && xs[i] <= 2) crossings.push_back((double(i) - (xs[i] - 2) / (xs[i] - xs[i - 1])) / 240.0);
    CHECK(crossings.size() >= 4);
    if (crossings.size() >= 4)
        CHECK_NEAR((crossings.back() - crossings.front()) / double(crossings.size() - 1), 0.5, 0.025);
}

TEST(damping_ratio_controls_how_quickly_the_oscillation_dies) {
    auto amplitude_after_1s = [](Real zeta) {
        const auto xs = spring_run(2.0f, zeta, 240);
        double peak = 0;
        for (size_t i = 120; i < xs.size(); ++i) peak = std::max(peak, std::fabs(xs[i] - 2));
        return peak;
    };
    const double none = amplitude_after_1s(0.0f), some = amplitude_after_1s(0.2f), lots = amplitude_after_1s(0.7f);
    CHECK(none > 0.25);        // an undamped spring keeps nearly its full 0.3 m swing
    CHECK(some < none * 0.7);
    CHECK(lots < some);

    // Critical damping: returns without overshooting.
    const auto critical = spring_run(2.0f, 1.0f, 240 * 2);
    double lowest = 10;
    for (double x : critical) lowest = std::min(lowest, x);
    CHECK(lowest > 2 - 0.01);
    CHECK_NEAR(critical.back(), 2, 0.005);
}

// ---- prismatic ------------------------------------------------------------------------------------

namespace {

// A crate on a rail tilted 30 degrees downhill (axis points down the slope).
struct Rail {
    World w;
    Vec2 axis{std::cos(-0.5236f), std::sin(-0.5236f)};
    Vec2 origin{0, 6};
    int body = 0;
    explicit Rail(bool limit = false, Real lower = 0, Real upper = 0) {
        body = w.add(plank(0.3f, 0.3f, origin));
        Joint j = Joint::prismatic(w.bodies, -1, body, origin, axis);
        j.enable_limit = limit;
        j.lower = lower;
        j.upper = upper;
        w.add_joint(j);
    }
    double position() const { return double(dot(w.bodies[static_cast<size_t>(body)].pos - origin, axis)); }
    double off_rail() const { return double(dot(perp(axis), w.bodies[static_cast<size_t>(body)].pos - origin)); }
};

}  // namespace

TEST(box_slides_down_a_rail_at_g_sin_theta) {
    Rail rail;
    run(rail.w, 180, 1.0f / 120.0f);  // 1.5 s
    const double expected = 0.5 * kG * 0.5 * 1.5 * 1.5;  // a = g sin(30) along the rail
    CHECK_NEAR(rail.position(), expected, 0.02 * expected);
    CHECK_NEAR(rail.off_rail(), 0, 0.01);               // it stays on the rail
    CHECK_NEAR(rail.w.bodies[0].angle, 0, 0.01);        // and does not twist
}

TEST(prismatic_limits_stop_the_box_at_both_ends) {
    Rail down(true, -1.0f, 2.0f);
    run(down.w, 360, 1.0f / 120.0f);
    CHECK_NEAR(down.position(), 2.0, 0.04);                     // gravity holds it against the upper limit
    CHECK_NEAR(double(down.w.bodies[0].vel.length()), 0, 0.05);

    Rail up(true, -1.0f, 2.0f);
    up.w.bodies[0].vel = up.axis * -8.0f;  // fired up the rail
    double lowest = 10;
    for (int i = 0; i < 240; ++i) {
        up.w.step(1.0f / 120.0f);
        lowest = std::min(lowest, up.position());
    }
    CHECK(lowest < -0.9);   // it reached the lower limit
    CHECK(lowest > -1.06);  // and was stopped there
}

TEST(prismatic_motor_drives_the_box_at_its_speed) {
    World w;  // horizontal rail; gravity must not pull the box off it
    w.add(plank(0.3f, 0.3f, {0, 3}));
    Joint j = Joint::prismatic(w.bodies, -1, 0, {0, 3}, {1, 0});
    j.enable_motor = true;
    j.motor_speed = 2;
    j.max_motor = 1000;
    w.add_joint(j);
    run(w, 120, 1.0f / 120.0f);
    CHECK_NEAR(w.bodies[0].vel.x, 2, 0.02);
    CHECK_NEAR(w.bodies[0].pos.y, 3, 0.01);
    CHECK_NEAR(w.bodies[0].vel.y, 0, 0.01);
}

// ---- mouse ----------------------------------------------------------------------------------------

TEST(mouse_joint_drags_a_body_to_the_target_without_overshoot) {
    World w;
    w.gravity = {};
    w.add(plank(0.3f, 0.3f, {0, 0}));
    Joint j = Joint::mouse(w.bodies, 0, {0, 0});
    j.target = {3, 2};
    j.damping_ratio = 1;
    w.add_joint(j);

    double overshoot = 0;
    for (int i = 0; i < 480; ++i) {
        w.step(1.0f / 120.0f);
        overshoot = std::max(overshoot, double(w.bodies[0].pos.x) - 3.0);
    }
    CHECK_NEAR(w.bodies[0].pos.x, 3, 0.02);
    CHECK_NEAR(w.bodies[0].pos.y, 2, 0.02);
    CHECK_NEAR(w.bodies[0].vel.length(), 0, 0.02);
    CHECK(overshoot < 0.05);
}

TEST(mouse_joint_force_is_capped) {
    World w;
    w.gravity = {};
    w.add(plank(0.3f, 0.3f, {0, 0}));
    const double mass = double(w.bodies[0].mass);
    Joint j = Joint::mouse(w.bodies, 0, {0, 0});
    j.target = {50, 0};  // very far: an uncapped spring would hurl the body
    j.max_force = static_cast<Real>(5.0 * mass);  // 5 m/s^2
    w.add_joint(j);
    run(w, 24, 1.0f / 120.0f);  // 0.2 s
    CHECK(double(w.bodies[0].vel.x) <= 5.0 * 0.2 * 1.05);
    CHECK(double(w.bodies[0].vel.x) > 5.0 * 0.2 * 0.8);  // but it is being pulled at about the cap
}

// ---- joints inside the World ----------------------------------------------------------------------

TEST(jointed_bodies_do_not_collide_unless_asked) {
    auto contacts_after_step = [](bool collide_connected) {
        World w;
        w.gravity = {};
        w.add(plank(0.5f, 0.5f, {0, 0}));
        w.add(plank(0.5f, 0.5f, {0.8f, 0}));  // overlapping by 0.2
        Joint j = Joint::revolute(w.bodies, 0, 1, {0.4f, 0});
        j.collide_connected = collide_connected;
        w.add_joint(j);
        w.step(1.0f / 60.0f);
        return w.contacts().size();
    };
    CHECK(contacts_after_step(false) == 0);
    CHECK(contacts_after_step(true) == 1);
}

TEST(truncate_removes_only_the_joints_that_lose_a_body) {
    World w;
    for (int i = 0; i < 4; ++i) w.add(disc(0.2f, {static_cast<Real>(i), 0}));
    w.add_joint(Joint::distance(w.bodies, 0, 1, {0, 0}, {1, 0}));
    w.add_joint(Joint::distance(w.bodies, 2, 3, {2, 0}, {3, 0}));
    w.add_joint(Joint::revolute(w.bodies, -1, 1, {1, 0}));  // to the world: keeps working
    w.truncate(3);
    CHECK(w.joints.size() == 2);  // the one touching body 3 is gone
    w.step(1.0f / 60.0f);          // and nothing dangles
    CHECK(w.joints[0].b == 1 && w.joints[1].a == -1);
}

// A rope pendulum swings down and strikes a crate on a ledge. The rope (a joint) and the impact (a
// contact) are solved in the same sweeps: the rope must hold while the crate is knocked away.
TEST(rope_pendulum_knocks_a_crate_off_while_the_rope_holds) {
    World w;
    Body ledge(Shape::make_polygon(Polygon::box(3, 0.1f)), {0, 1.5f}, 0, BodyType::Static);
    w.add(ledge);
    w.add(plank(0.4f, 0.4f, {0, 2.0f}));      // crate on the ledge, its left face at x = -0.4
    const Vec2 pivot{0, 5};
    w.add(disc(0.3f, {-3, 5}));               // ball level with the pivot, 3 m of rope out
    w.bodies[1].restitution = 0;
    w.bodies[2].restitution = 0;
    w.add_joint(Joint::distance(w.bodies, -1, 2, pivot, {-3, 5}));

    double worst_rope = 0;
    for (int i = 0; i < 360; ++i) {
        w.step(1.0f / 120.0f);
        worst_rope = std::max(worst_rope, std::fabs(double(distance(pivot, w.bodies[2].pos)) - 3.0));
    }
    CHECK(worst_rope < 0.03);
    CHECK(w.bodies[1].pos.x > 0.3);  // the crate was struck and shoved to the right
}

// A joint that starts out violated must pull itself back: a body knocked off its rail, and twisted,
// returns to the rail and to its reference angle.
TEST(prismatic_joint_corrects_being_off_the_rail_and_twisted) {
    World w;
    w.gravity = {};
    w.add(plank(0.3f, 0.3f, {0, 3}));
    w.add_joint(Joint::prismatic(w.bodies, -1, 0, {0, 3}, {1, 0}));
    w.bodies[0].pos.y = 3.25f;     // 25 cm off the rail
    w.bodies[0].set_angle(0.25f);  // and twisted by 0.25 rad
    run(w, 240, 1.0f / 120.0f);
    CHECK_NEAR(w.bodies[0].pos.y, 3, 0.01);
    CHECK_NEAR(w.bodies[0].angle, 0, 0.01);
}

// Same for the hinge and the rod: displaced anchors are drawn back together.
TEST(revolute_and_distance_joints_correct_an_initial_error) {
    World w;
    w.gravity = {};
    w.add(disc(0.3f, {0, 0}));
    w.add(disc(0.3f, {2, 0}));
    w.add_joint(Joint::revolute(w.bodies, -1, 0, {0, 0}));
    w.add_joint(Joint::distance(w.bodies, 0, 1, {0, 0}, {2, 0}));
    w.bodies[0].pos = {0.2f, -0.1f};  // hinge pulled apart
    w.bodies[1].pos = {2.5f, 0};      // rod stretched
    run(w, 360, 1.0f / 120.0f);
    CHECK(anchor_error(w, w.joints[0]) < 0.01);                                  // hinge: anchors coincide
    CHECK_NEAR(double(distance(w.bodies[0].pos, w.bodies[1].pos)), 2.0, 0.01);  // rod: back to 2 m
}

// ---- heavy loads and error correction: the reason joints have their own position pass -------------

namespace {

// An 8-link chain (links 0.05 kg) hanging from the world with a ball of the given density on the end.
// Returns the worst gap between any pair of hinge anchors over `seconds`.
double worst_hinge_gap(Real ball_density, int sweeps, bool warm, Real ball_speed, double seconds) {
    World w;
    w.solver.iterations = sweeps;
    w.solver.warm_starting = warm;
    const Vec2 top{0, 8};
    const int links = 8;
    for (int i = 0; i < links; ++i) w.add(plank(0.05f, 0.25f, {top.x, top.y - 0.25f - 0.5f * static_cast<Real>(i)}));
    const int ball = w.add(Body(Shape::make_circle(0.5f), {top.x, top.y - 0.5f * links - 0.5f}, 0, BodyType::Dynamic, ball_density));
    w.add_joint(Joint::revolute(w.bodies, -1, 0, top));
    for (int i = 1; i < links; ++i) w.add_joint(Joint::revolute(w.bodies, i - 1, i, {top.x, top.y - 0.5f * static_cast<Real>(i)}));
    w.add_joint(Joint::revolute(w.bodies, links - 1, ball, {top.x, top.y - 0.5f * links}));
    w.bodies[static_cast<size_t>(ball)].vel = {ball_speed, 0};

    double worst = 0;
    for (int i = 0; i < static_cast<int>(seconds * 120); ++i) {
        w.step(1.0f / 120.0f);
        for (const Joint& j : w.joints) worst = std::max(worst, anchor_error(w, j));
        if (worst > 50) break;  // torn apart; no need to keep going
    }
    return worst;
}

}  // namespace

// Weight hanging from a light chain, at rest. The load is 16x to 314x heavier than a link. Warm starting
// carries the tension from step to step, so the chain simply holds; without it ten sweeps cannot push
// the weight's tension all the way up the chain in one step and the links stretch apart.
TEST(a_heavy_weight_hangs_from_a_light_chain) {
    CHECK(worst_hinge_gap(1.0f, 10, true, 0, 6) < 0.01);    // ball 16x a link
    CHECK(worst_hinge_gap(20.0f, 10, true, 0, 6) < 0.1);    // ball 314x a link
    const double cold = worst_hinge_gap(5.0f, 10, false, 0, 6);  // ball 79x a link, no warm start
    const double warm = worst_hinge_gap(5.0f, 10, true, 0, 6);
    CHECK(warm < 0.05);
    CHECK(cold > 20 * warm);
}

// Swung hard, the same 314x ball stays attached. (With only 4 sweeps it would not: each sweep moves
// information one joint along the chain.)
TEST(a_swung_heavy_weight_is_held) {
    CHECK(worst_hinge_gap(20.0f, 10, true, 3, 6) < 0.2);
}

namespace {

// The rope bridge from the demo: 12 planks laid along a sagging curve, both ends pinned to the world.
struct Bridge {
    World w;
    std::vector<int> link;
    Bridge() {
        constexpr int kLinks = 12;
        const Vec2 left{1.0f, 6.4f}, right{7.4f, 6.4f};
        std::vector<Vec2> pts;
        for (int i = 0; i <= kLinks; ++i) {
            const Real t = static_cast<Real>(i) / kLinks;
            pts.push_back({left.x + (right.x - left.x) * t, left.y - 1.4f * 4 * t * (1 - t)});
        }
        for (int i = 0; i < kLinks; ++i) {
            const Vec2 d = pts[static_cast<size_t>(i) + 1] - pts[static_cast<size_t>(i)];
            Body b(Shape::make_polygon(Polygon::box(d.length() * 0.5f, 0.07f)),
                   (pts[static_cast<size_t>(i)] + pts[static_cast<size_t>(i) + 1]) * 0.5f, std::atan2(d.y, d.x));
            b.restitution = 0;
            link.push_back(w.add(b));
        }
        w.add_joint(Joint::revolute(w.bodies, -1, link.front(), pts.front()));
        for (int i = 1; i < kLinks; ++i)
            w.add_joint(Joint::revolute(w.bodies, link[static_cast<size_t>(i) - 1], link[static_cast<size_t>(i)],
                                        pts[static_cast<size_t>(i)]));
        w.add_joint(Joint::revolute(w.bodies, link.back(), -1, pts.back()));
    }
    double worst_gap() const {
        double g = 0;
        for (const Joint& j : w.joints)
            if (j.type == JointType::Revolute) g = std::max(g, anchor_error(w, j));  // hinges only: a mouse joint is MEANT to stretch
        return g;
    }
};

}  // namespace

// Both ends pinned, sagging under its own weight, then the middle yanked 2.2 m down by a mouse joint.
TEST(a_pinned_rope_bridge_holds_together_while_dragged) {
    Bridge bridge;
    double worst = 0;
    for (int i = 0; i < 300; ++i) {
        bridge.w.step(1.0f / 120.0f);
        worst = std::max(worst, bridge.worst_gap());
    }
    CHECK(worst < 0.01);  // settling under gravity alone
    const int mid = bridge.link[5];
    Joint grab = Joint::mouse(bridge.w.bodies, mid, bridge.w.bodies[static_cast<size_t>(mid)].pos);
    grab.target = bridge.w.bodies[static_cast<size_t>(mid)].pos + Vec2{0, -2.2f};
    bridge.w.add_joint(grab);
    for (int i = 0; i < 240; ++i) {
        bridge.w.step(1.0f / 120.0f);
        worst = std::max(worst, bridge.worst_gap());
    }
    CHECK(worst < 0.02);
}

// Position correction must never make a chain worse. In a chain every link belongs to two joints; if each
// joint computed its fix from the same starting poses, the shared link would receive both fixes at once
// and overshoot, turning a small error into a bigger one.
TEST(position_correction_reduces_error_on_a_chain) {
    Bridge bridge;
    bridge.w.bodies[static_cast<size_t>(bridge.link[6])].pos += Vec2{0.05f, -0.04f};  // knock one plank out of place
    const double before = bridge.worst_gap();
    SolverSettings settings;
    settings.joint_position_iterations = 1;
    solve_joint_positions(bridge.w.bodies, bridge.w.joints, settings);
    const double after_one = bridge.worst_gap();
    CHECK(after_one < before);
    for (int i = 0; i < 20; ++i) solve_joint_positions(bridge.w.bodies, bridge.w.joints, settings);
    CHECK(bridge.worst_gap() < 0.25 * before);  // and keeps shrinking
}

// A badly violated joint is eased back, not teleported: one step removes at most 20 cm per correction
// pass, so a 1.5 m error is still mostly there after one step but is gone after a while.
TEST(large_joint_errors_are_eased_back_not_teleported) {
    World w;
    w.gravity = {};
    w.add(disc(0.3f, {0, 0}));
    w.add_joint(Joint::revolute(w.bodies, -1, 0, {0, 0}));
    w.bodies[0].pos = {1.5f, 0};
    w.step(1.0f / 120.0f);
    CHECK(anchor_error(w, w.joints[0]) > 1.5 - 4 * 0.2 - 0.01);  // at most 4 passes x 0.2 m
    CHECK(anchor_error(w, w.joints[0]) < 1.5);                    // but it did start coming back
    run(w, 240, 1.0f / 120.0f);
    CHECK(anchor_error(w, w.joints[0]) < 0.005);
}

// Correcting a joint's error moves bodies, it does not set them moving. A hinge pulled 15 cm off its pin
// in zero gravity comes back without ever picking up speed (feeding the error back as velocity, the old
// Baumgarte way, would fling the body at about 3.6 m/s).
TEST(joint_error_is_removed_without_adding_energy) {
    World w;
    w.gravity = {};
    w.add(disc(0.3f, {0, 0}));
    w.add_joint(Joint::revolute(w.bodies, -1, 0, {0, 0}));
    w.bodies[0].pos = {0.15f, 0};
    double fastest = 0;
    for (int i = 0; i < 60; ++i) {
        w.step(1.0f / 120.0f);
        fastest = std::max(fastest, double(w.bodies[0].vel.length()));
    }
    CHECK(fastest < 0.01);
    CHECK(anchor_error(w, w.joints[0]) < 0.005);
}

// A limit that is not reached must not reach out and pull the body toward it.
TEST(an_unreached_limit_does_not_pull_the_body) {
    World w;
    w.gravity = {};
    w.add(plank(0.3f, 0.3f, {0, 3}));
    Joint j = Joint::revolute(w.bodies, -1, 0, {0, 3});
    j.enable_limit = true;
    j.lower = -1.0f;
    j.upper = 1.0f;
    w.add_joint(j);
    w.bodies[0].set_angle(0.4f);  // comfortably inside [-1, 1]
    run(w, 240, 1.0f / 120.0f);
    CHECK_NEAR(w.bodies[0].angle, 0.4, 1e-3);
}
