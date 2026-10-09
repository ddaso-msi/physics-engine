#include "test.hpp"
#include <phys/body.hpp>

#include <vector>

using namespace phys;

namespace {
constexpr double kEps = 1e-4;
constexpr double kPiD = 3.14159265358979323846;

Polygon regular_polygon(int n, Real radius) {
    std::vector<Vec2> pts;
    for (int k = 0; k < n; ++k) {
        Real a = 2 * kPi * static_cast<Real>(k) / static_cast<Real>(n);
        pts.push_back({radius * std::cos(a), radius * std::sin(a)});
    }
    return *Polygon::from_points(pts);
}
}  // namespace

TEST(circle_mass) {
    MassData md = compute_mass(Shape::make_circle(2.0f), 3.0f);
    CHECK_NEAR(md.mass, 3.0 * kPiD * 4.0, 1e-3);
    CHECK_NEAR(md.inertia, 0.5 * md.mass * 4.0, 1e-3);
}

TEST(box_mass_matches_rectangle_formula) {
    // w = 4, h = 2, rho = 1.5: m = rho*w*h, I = m*(w^2 + h^2)/12
    MassData md = compute_mass(Shape::make_polygon(Polygon::box(2.0f, 1.0f)), 1.5f);
    CHECK_NEAR(md.mass, 12.0, kEps);
    CHECK_NEAR(md.inertia, 12.0 * (16.0 + 4.0) / 12.0, kEps);
}

// Regular n-gon with circumradius R: I = m R^2 / 6 * (1 + 2 cos^2(pi/n)).
TEST(regular_polygon_inertia) {
    for (int n : {3, 4, 6, 8}) {
        const double R = 1.5, rho = 2.0;
        MassData md = compute_mass(Shape::make_polygon(regular_polygon(n, 1.5f)), 2.0f);
        double area = 0.5 * n * R * R * std::sin(2 * kPiD / n);
        double c = std::cos(kPiD / n);
        CHECK_NEAR(md.mass, rho * area, 1e-3);
        CHECK_NEAR(md.inertia, md.mass * R * R / 6.0 * (1 + 2 * c * c), 1e-3);
    }
}

// A box defined far from the origin is recentred; its mass properties must not change.
TEST(polygon_recentres_on_centroid) {
    const Vec2 pts[4] = {{10, 10}, {12, 10}, {12, 13}, {10, 13}};  // 2 x 3, centroid (11, 11.5)
    auto poly = Polygon::from_points(pts);
    CHECK(poly.has_value());
    CHECK_NEAR(poly->centroid_shift.x, 11, kEps);
    CHECK_NEAR(poly->centroid_shift.y, 11.5, kEps);

    Vec2 sum;
    for (int i = 0; i < poly->count; ++i) sum += poly->vertices[i];
    CHECK_NEAR(sum.x, 0, kEps);
    CHECK_NEAR(sum.y, 0, kEps);

    MassData md = compute_mass(Shape::make_polygon(*poly), 1.0f);
    CHECK_NEAR(md.mass, 6.0, kEps);
    CHECK_NEAR(md.inertia, 6.0 * (4.0 + 9.0) / 12.0, kEps);
}

TEST(polygon_validation) {
    const Vec2 square_cw[4] = {{-1, -1}, {-1, 1}, {1, 1}, {1, -1}};
    auto cw = Polygon::from_points(square_cw);
    CHECK(cw.has_value());  // clockwise input is flipped
    if (cw) {
        MassData md = compute_mass(Shape::make_polygon(*cw), 1.0f);
        CHECK_NEAR(md.mass, 4.0, kEps);  // positive: winding was normalised
    }

    const Vec2 two[2] = {{0, 0}, {1, 0}};
    CHECK(!Polygon::from_points(two).has_value());

    const Vec2 arrow[4] = {{0, 0}, {2, 1}, {0, 2}, {1, 1}};  // concave
    CHECK(!Polygon::from_points(arrow).has_value());

    const Vec2 collinear[4] = {{0, 0}, {1, 0}, {2, 0}, {1, 1}};  // (0,0),(1,0),(2,0) in a line
    CHECK(!Polygon::from_points(collinear).has_value());

    std::vector<Vec2> nine;
    for (int k = 0; k < 9; ++k) {
        Real a = 2 * kPi * static_cast<Real>(k) / 9;
        nine.push_back({std::cos(a), std::sin(a)});
    }
    CHECK(!Polygon::from_points(nine).has_value());
}

TEST(polygon_normals_point_outward) {
    Polygon box = Polygon::box(2.0f, 1.0f);
    for (int i = 0; i < box.count; ++i) {
        CHECK_NEAR(box.normals[i].length(), 1, kEps);
        CHECK(dot(box.normals[i], box.vertices[i]) > 0);  // origin is inside, so outward means away from it
    }
    CHECK_NEAR(box.normals[0].y, -1, kEps);  // bottom edge
    CHECK_NEAR(box.normals[1].x, 1, kEps);   // right edge
}

TEST(static_body_has_infinite_mass) {
    Body ground(Shape::make_polygon(Polygon::box(10, 1)), {0, 0}, 0, BodyType::Static);
    CHECK(ground.inv_mass == 0 && ground.inv_inertia == 0);
    ground.apply_force({100, 100});
    ground.apply_impulse_at({50, 50}, {1, 1});
    ground.integrate(1.0f, {0, -9.81f});
    CHECK(ground.pos.x == 0 && ground.pos.y == 0);
    CHECK(ground.vel.length_sq() == 0 && ground.w == 0);
}

TEST(gravity_is_mass_independent) {
    Body light(Shape::make_circle(0.5f), {0, 0}, 0, BodyType::Dynamic, 1.0f);
    Body heavy(Shape::make_circle(0.5f), {0, 0}, 0, BodyType::Dynamic, 100.0f);
    for (int i = 0; i < 60; ++i) {
        light.integrate(1.0f / 60.0f, {0, -9.81f});
        heavy.integrate(1.0f / 60.0f, {0, -9.81f});
    }
    CHECK_NEAR(light.pos.y, heavy.pos.y, 1e-6);
    CHECK_NEAR(light.vel.y, -9.81, 1e-3);
}

// A force applied at an offset point accelerates the centre of mass by F/m (wherever it acts) and
// spins the body with torque r x F. Constant torque: w = (tau / I) t exactly for this integrator.
TEST(force_at_point_gives_linear_and_angular_acceleration) {
    Body b(Shape::make_circle(1.0f), {0, 0}, 0);  // m = pi, I = pi/2
    const int n = 1000;
    const Real dt = 1.0f / static_cast<Real>(n);
    for (int i = 0; i < n; ++i) {
        b.apply_force_at({0, 10}, b.pos + Vec2{2, 0});  // r = (2,0) from the centre, so tau = 20
        b.integrate(dt, {});
    }
    CHECK_NEAR(b.vel.y, 10.0 / kPiD, 1e-2);
    CHECK_NEAR(b.vel.x, 0, kEps);
    CHECK_NEAR(b.w, 20.0 / (kPiD / 2), 1e-2);
}

TEST(forces_are_cleared_after_each_step) {
    Body b(Shape::make_circle(1.0f), {0, 0}, 0);
    b.apply_force({10, 0});
    b.apply_torque(5);
    b.integrate(0.1f, {});
    Vec2 v = b.vel;
    Real w = b.w;
    b.integrate(0.1f, {});  // nothing applied this step: velocities unchanged
    CHECK_NEAR(b.vel.x, v.x, 1e-7);
    CHECK_NEAR(b.w, w, 1e-7);
}

// Angular momentum about the world origin, L = m (p x v) + I w, changes by exactly p_hit x j.
TEST(impulse_changes_angular_momentum_by_r_cross_j) {
    Body b(Shape::make_polygon(Polygon::box(1.0f, 0.5f)), {3, 1}, 0.3f);
    auto L = [&] { return b.mass * cross(b.pos, b.vel) + b.inertia * b.w; };
    const Vec2 hit{4, 1.5f}, j{0, 2};
    Real before = L();
    b.apply_impulse_at(j, hit);
    CHECK_NEAR(L() - before, cross(hit, j), 1e-4);
    CHECK_NEAR(b.vel.y, 2 * b.inv_mass, 1e-5);
}

TEST(velocity_at_point_adds_rotation) {
    Body b(Shape::make_circle(1.0f), {1, 1}, 0);
    b.vel = {2, 0};
    b.w = 3;
    Vec2 v = b.velocity_at({2, 1});  // 1 m to the right of the centre: w x r = (0, 3)
    CHECK_NEAR(v.x, 2, kEps);
    CHECK_NEAR(v.y, 3, kEps);
}

TEST(free_spin_conserves_energy_and_advances_angle) {
    Body b(Shape::make_polygon(Polygon::box(1.0f, 0.5f)), {0, 0}, 0);
    b.w = 3;
    const Real ke = b.kinetic_energy();
    const int n = 1000;
    for (int i = 0; i < n; ++i) b.integrate(1.0f / static_cast<Real>(n), {});
    CHECK_NEAR(b.kinetic_energy(), ke, 1e-4);
    CHECK_NEAR(b.w, 3, 1e-6);
    CHECK_NEAR(b.q.c, std::cos(3.0), 1e-3);  // compared via cos/sin: the stored angle wraps at pi
    CHECK_NEAR(b.q.s, std::sin(3.0), 1e-3);
    CHECK(std::fabs(b.angle) <= kPi + 1e-5f);
}

TEST(contains_respects_rotation) {
    Body b(Shape::make_polygon(Polygon::box(1.0f, 0.5f)), {5, 5}, kPi / 2);  // now 1 wide, 2 tall
    CHECK(b.contains({5, 5.9f}));
    CHECK(!b.contains({5.9f, 5}));
    CHECK(!b.contains({7, 7}));

    Body c(Shape::make_circle(1.0f), {0, 0}, 0);
    CHECK(c.contains({0.7f, 0.7f}));
    CHECK(!c.contains({0.8f, 0.8f}));
}
