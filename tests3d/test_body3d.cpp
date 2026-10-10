#include "test.hpp"
#include <phys3d/body.hpp>

#include <algorithm>
#include <cstdint>

using namespace phys3d;

namespace {

constexpr double kG = 9.81;
constexpr double kPiD = 3.14159265358979323846;

struct Lcg {
    std::uint32_t s = 9;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Vec3 vec(Real extent = 3) { return {range(-extent, extent), range(-extent, extent), range(-extent, extent)}; }
    Quat rotation() {
        Quat q;
        do q = Quat{range(-1, 1), range(-1, 1), range(-1, 1), range(-1, 1)};
        while (q.length() < 0.2f);
        return q.normalized();
    }
};

#define CHECK_VEC(a, b, eps)                         \
    do {                                             \
        const Vec3 va_ = (a), vb_ = (b);             \
        CHECK_NEAR(va_.x, vb_.x, eps);               \
        CHECK_NEAR(va_.y, vb_.y, eps);               \
        CHECK_NEAR(va_.z, vb_.z, eps);               \
    } while (0)

double angle_between(Vec3 a, Vec3 b) {
    const double c = double(dot(a, b)) / (double(a.length()) * double(b.length()));
    return std::acos(std::max(-1.0, std::min(1.0, c)));
}

// An asymmetric box: three different principal moments (0.128, 0.333, 0.435).
Body brick() { return Body::solid_box({1.0f, 0.6f, 0.2f}, 1, {0, 0, 0}); }

// What a torque-free run does to the conserved quantities.
struct Drift {
    double energy_min = 1e9, energy_max = 0, energy_final = 0;  // as a fraction of the start
    double momentum_angle = 0;                                    // worst turn of L, radians
    double momentum_size = 0;                                     // worst relative change of |L|
};
Drift tumble(Body b, Gyroscopic mode, Real dt, double seconds) {
    const Vec3 l0 = b.angular_momentum();
    const double e0 = double(b.kinetic_energy());
    Drift d;
    for (int i = 0; i < static_cast<int>(seconds / double(dt)); ++i) {
        b.integrate(dt, {}, mode);
        const double e = double(b.kinetic_energy()) / e0;
        d.energy_min = std::min(d.energy_min, e);
        d.energy_max = std::max(d.energy_max, e);
        const Vec3 l = b.angular_momentum();
        d.momentum_angle = std::max(d.momentum_angle, angle_between(l, l0));
        d.momentum_size = std::max(d.momentum_size, std::fabs(double(l.length()) / double(l0.length()) - 1.0));
    }
    d.energy_final = double(b.kinetic_energy()) / e0;
    return d;
}

}  // namespace

// ---- mass properties and the immovable body -------------------------------------------------------

TEST(solid_shapes_get_their_mass_and_inertia) {
    const Body box = Body::solid_box({1, 2, 3}, 0.5f, {0, 0, 0});  // 2 x 4 x 6 = 48 volume
    CHECK_NEAR(box.mass, 24, 1e-3);
    CHECK_NEAR(box.inv_mass, 1.0 / 24, 1e-6);
    CHECK_NEAR(box.inertia_body.at(0, 0), 24.0 * (4 + 9) / 3, 1e-2);  // m (b^2 + c^2) / 3
    CHECK_NEAR(box.inv_inertia_body.at(0, 0) * box.inertia_body.at(0, 0), 1, 1e-4);

    const Body ball = Body::solid_sphere(2, 3, {0, 0, 0});
    CHECK_NEAR(ball.mass, 3.0 * 4.0 / 3.0 * kPiD * 8, 1e-2);
    CHECK_NEAR(ball.inertia_body.at(1, 1), 0.4 * double(ball.mass) * 4, 1e-2);
}

TEST(a_fixed_body_cannot_be_moved_or_spun) {
    Body wall = Body::fixed({1, 2, 3}, Quat::from_axis_angle({0, 1, 0}, 0.4f));
    CHECK(wall.inv_mass == 0 && wall.inv_inertia_world().cx.length_sq() == 0);
    wall.apply_force_at({100, 50, -20}, {4, 4, 4});
    wall.apply_impulse_at({30, 0, 0}, {0, 5, 0});
    wall.integrate(1.0f, {0, static_cast<Real>(-kG), 0});
    CHECK_VEC(wall.pos, Vec3(1, 2, 3), 0);
    CHECK(wall.vel.length_sq() == 0 && wall.w.length_sq() == 0);
    CHECK(wall.force.length_sq() == 0 && wall.torque.length_sq() == 0);  // accumulators still cleared

    // Even a velocity written into it by hand does not move it.
    wall.vel = {4, 0, 0};
    wall.w = {0, 3, 0};
    const Quat before = wall.q;
    wall.integrate(1.0f, {});
    CHECK_VEC(wall.pos, Vec3(1, 2, 3), 0);
    CHECK(wall.q.x == before.x && wall.q.y == before.y && wall.q.z == before.z && wall.q.w == before.w);
}

TEST(world_inverse_inertia_inverts_world_inertia_in_any_orientation) {
    Lcg rng;
    for (int i = 0; i < 100; ++i) {
        Body b = brick();
        b.q = rng.rotation();
        const Mat3 product = b.inertia_world() * b.inv_inertia_world();
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) CHECK_NEAR(product.at(r, c), r == c ? 1 : 0, 1e-3);
    }
}

// ---- linear motion and forces ---------------------------------------------------------------------

TEST(free_fall_in_3d) {
    Body light = Body::solid_sphere(0.2f, 1, {1, 10, -2}), heavy = Body::solid_sphere(0.2f, 500, {1, 10, -2});
    const int n = 120;
    const Real dt = 1.0f / 120.0f;
    for (int i = 0; i < n; ++i) {
        light.integrate(dt, {0, static_cast<Real>(-kG), 0});
        heavy.integrate(dt, {0, static_cast<Real>(-kG), 0});
    }
    // Semi-implicit Euler falls g dt^2 n(n+1)/2, a hair more than g t^2 / 2 (see the 2D stage 2 tests).
    CHECK_NEAR(light.pos.y, 10 - kG * double(dt) * double(dt) * n * (n + 1) / 2, 1e-3);
    CHECK_NEAR(light.vel.y, -kG, 1e-3);
    CHECK_NEAR(light.pos.y, heavy.pos.y, 1e-5);  // gravity is mass independent
    CHECK_NEAR(light.pos.x, 1, 1e-6);
    CHECK_NEAR(light.pos.z, -2, 1e-6);
}

TEST(force_at_a_point_pushes_and_twists) {
    // A force along +y applied 2 m out along +x: torque r x F = (2,0,0) x (0,10,0) = (0,0,20), about z.
    Body b = brick();
    const double izz = double(b.inertia_body.at(2, 2));
    const int n = 1000;
    const Real dt = 1.0f / static_cast<Real>(n);
    for (int i = 0; i < n; ++i) {
        b.apply_force_at({0, 10, 0}, b.pos + Vec3{2, 0, 0});
        b.integrate_velocity(dt, {});  // velocity only: keep the body where it is so r stays (2,0,0)
    }
    CHECK_VEC(b.vel, Vec3(0, static_cast<Real>(10.0 / double(b.mass)), 0), 1e-2);  // F/m wherever it acts
    CHECK_VEC(b.w, Vec3(0, 0, static_cast<Real>(20.0 / izz)), 0.05);               // tau / I about a principal axis
}

TEST(forces_and_torques_are_cleared_after_each_step) {
    Body b = brick();
    b.apply_force({5, 0, 0});
    b.apply_torque({0, 0, 3});
    b.integrate(0.1f, {});
    const Vec3 v = b.vel, w = b.w;
    b.integrate(0.1f, {}, Gyroscopic::Off);  // nothing applied: with the gyroscopic term off, nothing changes
    CHECK_VEC(b.vel, v, 1e-7);
    CHECK_VEC(b.w, w, 1e-7);
}

// Linear momentum changes by j, and angular momentum about the world origin, L = m (p x v) + I w,
// changes by exactly (hit point) x j, whatever the body's orientation.
TEST(impulse_changes_momentum_and_angular_momentum_exactly) {
    Lcg rng;
    for (int i = 0; i < 100; ++i) {
        Body b = brick();
        b.pos = rng.vec(4);
        b.q = rng.rotation();
        b.vel = rng.vec(2);
        b.w = rng.vec(2);
        auto momentum = [&] { return b.vel * b.mass; };
        auto about_origin = [&] { return cross(b.pos, b.vel) * b.mass + b.angular_momentum(); };
        const Vec3 hit = b.pos + rng.vec(1), j = rng.vec(2);
        const Vec3 p0 = momentum(), l0 = about_origin();
        b.apply_impulse_at(j, hit);
        CHECK_VEC(momentum() - p0, j, 1e-3);
        CHECK_VEC(about_origin() - l0, cross(hit, j), 5e-3);
    }
}

TEST(velocity_at_a_point_adds_the_spin) {
    Body b = brick();
    b.pos = {1, 1, 1};
    b.vel = {2, 0, 0};
    b.w = {0, 0, 3};  // spinning about z
    CHECK_VEC(b.velocity_at({2, 1, 1}), Vec3(2, 3, 0), 1e-5);  // 1 m out along x: w x r = (0, 3, 0)
    CHECK_VEC(b.velocity_at({1, 1, 5}), Vec3(2, 0, 0), 1e-5);  // on the spin axis: no extra velocity
}

TEST(kinetic_energy_counts_translation_and_rotation) {
    Body b = brick();
    b.vel = {3, 0, 0};
    b.w = {0, 2, 0};
    CHECK_NEAR(b.kinetic_energy(), 0.5 * double(b.mass) * 9 + 0.5 * double(b.inertia_body.at(1, 1)) * 4, 1e-3);
}

// ---- rotation without torque: the gyroscopic term --------------------------------------------------

TEST(spin_about_a_principal_axis_is_steady_and_turns_the_body_at_that_rate) {
    for (Gyroscopic mode : {Gyroscopic::Off, Gyroscopic::Explicit, Gyroscopic::Implicit}) {
        Body b = brick();
        b.w = {0, 0, 2.5f};
        for (int i = 0; i < 240; ++i) b.integrate(1.0f / 240.0f, {}, mode);
        CHECK_VEC(b.w, Vec3(0, 0, 2.5f), 1e-4);
        const Quat expected = Quat::from_axis_angle({0, 0, 1}, 2.5f);  // 1 s at 2.5 rad/s
        CHECK_VEC(rotate(b.q, Vec3(1, 0, 0)), rotate(expected, Vec3(1, 0, 0)), 2e-3);
    }
}

TEST(a_sphere_needs_no_gyroscopic_term) {
    // Its inertia is the same about every axis, so I w is parallel to w and w x (I w) vanishes.
    for (Gyroscopic mode : {Gyroscopic::Off, Gyroscopic::Explicit, Gyroscopic::Implicit}) {
        Body b = Body::solid_sphere(0.5f, 1, {0, 0, 0});
        b.w = {3, -2, 1};
        for (int i = 0; i < 600; ++i) b.integrate(1.0f / 120.0f, {}, mode);
        CHECK_VEC(b.w, Vec3(3, -2, 1), 1e-3);
    }
}

// Without the gyroscopic term w never changes, so L = I_world w turns with the body: angular momentum
// appears from nowhere. With it, L holds its direction.
TEST(a_tumbling_box_conserves_angular_momentum_only_with_the_gyroscopic_term) {
    Body b = brick();
    b.w = {3, 2, 1};
    const Drift off = tumble(b, Gyroscopic::Off, 1.0f / 240.0f, 10);
    CHECK(off.momentum_angle > 0.5);  // L swung through more than half a radian

    const Drift implicit = tumble(b, Gyroscopic::Implicit, 1.0f / 240.0f, 10);
    CHECK(implicit.momentum_angle < 0.03);
    CHECK(implicit.momentum_size < 0.1);

    // The explicit form gets the direction right too (its fault is energy, tested below).
    CHECK(tumble(b, Gyroscopic::Explicit, 1.0f / 240.0f, 10).momentum_angle < 0.03);
}

// The inertia that resists a torque is the WORLD-frame one. Turn the brick a quarter turn about x, so its
// own y axis (moment 0.333) now points along world z, and twist it about world z: it must respond with
// the y moment, not the z moment (0.435) it would use if the tensor were not rotated with the body.
TEST(torque_on_a_turned_body_meets_the_world_frame_inertia) {
    Body b = Body::solid_box({1.0f, 0.6f, 0.2f}, 1, {0, 0, 0}, Quat::from_axis_angle({1, 0, 0}, kPi / 2));
    const double iyy = double(b.inertia_body.at(1, 1)), izz = double(b.inertia_body.at(2, 2));
    for (int i = 0; i < 1000; ++i) {
        b.apply_torque({0, 0, 2});
        b.integrate_velocity(1.0f / 1000.0f, {});  // velocity only, so the body does not turn away
    }
    CHECK_NEAR(b.w.z, 2.0 / iyy, 0.02);
    CHECK(std::fabs(2.0 / iyy - 2.0 / izz) > 1.0);  // the two candidate answers are far apart
    CHECK_NEAR(b.w.x, 0, 1e-3);
    CHECK_NEAR(b.w.y, 0, 1e-3);
}

// The explicit form pumps energy in; the implicit form can only lose a little. On a fast-spinning thin
// plate the explicit form explodes within seconds.
TEST(explicit_gyroscopics_gain_energy_and_implicit_never_does) {
    Body b = brick();
    b.w = {3, 2, 1};
    const Drift explicit_run = tumble(b, Gyroscopic::Explicit, 1.0f / 60.0f, 10);
    CHECK(explicit_run.energy_final > 1.3);
    CHECK(explicit_run.energy_min >= 1.0 - 1e-4);  // it only ever goes up

    const Drift implicit_run = tumble(b, Gyroscopic::Implicit, 1.0f / 60.0f, 10);
    CHECK(implicit_run.energy_max <= 1.0 + 1e-4);  // it only ever goes down...
    CHECK(implicit_run.energy_final > 0.6);         // ...and not catastrophically (29% in 10 s at this coarse step)
    const Drift implicit_fine = tumble(b, Gyroscopic::Implicit, 1.0f / 240.0f, 10);
    CHECK(implicit_fine.energy_final > 0.88);       // a finer step loses less (9%)
    CHECK(implicit_fine.energy_final > implicit_run.energy_final);

    Body plate = Body::solid_box({1, 1, 0.05f}, 1, {0, 0, 0});
    plate.w = {0.5f, 0.3f, 30};
    const Drift plate_explicit = tumble(plate, Gyroscopic::Explicit, 1.0f / 60.0f, 5);
    CHECK(!(plate_explicit.energy_max < 100));  // blown up (or NaN)
    const Drift plate_implicit = tumble(plate, Gyroscopic::Implicit, 1.0f / 60.0f, 5);
    CHECK(plate_implicit.energy_max <= 1.0 + 1e-4);
    CHECK(plate_implicit.energy_final > 0.99);
}

// A symmetric body (I1 = I2 != I3) spun mostly about its symmetry axis: Euler's equations have a closed
// form. Seen from the body, w circles the symmetry axis at the rate Omega = (I3 - I1) / I1 * w3.
TEST(symmetric_top_precesses_at_the_closed_form_rate) {
    Body top = Body::solid_box({0.5f, 0.5f, 1.0f}, 1, {0, 0, 0});
    const double i1 = double(top.inertia_body.at(0, 0)), i3 = double(top.inertia_body.at(2, 2));
    top.w = {1, 0, 4};
    const double omega = (i3 - i1) / i1 * 4.0, seconds = 2.0;
    const Real dt = 1.0f / 1000.0f;
    for (int i = 0; i < 2000; ++i) top.integrate(dt, {}, Gyroscopic::Implicit);
    const Vec3 wb = inv_rotate(top.q, top.w);
    CHECK_NEAR(wb.x, std::cos(omega * seconds), 0.02);
    CHECK_NEAR(wb.y, std::sin(omega * seconds), 0.02);
    CHECK_NEAR(wb.z, 4, 1e-3);  // the spin about the symmetry axis itself never changes
    CHECK(std::fabs(omega) > 1.0);  // the test really is about a moving target
}

// The tennis racket (Dzhanibekov) effect: rotation about the axis of the middle moment of inertia is
// unstable. A 1% wobble grows until the body flips right over; about the other two axes it stays put.
TEST(spin_about_the_intermediate_axis_flips_over) {
    auto lowest_spin_about = [](int axis) {
        Body b = brick();  // moments 0.128 < 0.333 < 0.435 about x, y, z
        Vec3 w{0.05f, 0.05f, 0.05f};
        w[axis] = 5;
        b.w = w;
        double lowest = 1e9;
        for (int i = 0; i < 30 * 240; ++i) {
            b.integrate(1.0f / 240.0f, {}, Gyroscopic::Implicit);
            lowest = std::min(lowest, double(inv_rotate(b.q, b.w)[axis]));
        }
        return lowest;
    };
    CHECK(lowest_spin_about(0) > 4.9);   // smallest moment: stable
    CHECK(lowest_spin_about(2) > 4.9);   // largest moment: stable
    CHECK(lowest_spin_about(1) < -4.0);  // the middle one: it turned completely over
}

// The gyroscopic term is a torque at right angles to w, so it can do no work. At extreme spin, several
// radians per step, the single Newton iteration of the implicit step overshoots, and without a guard it
// hands back far more rotational energy than it was given (a stress test saw x780 in one step).
TEST(gyroscopic_step_never_adds_energy_even_at_extreme_spin) {
    Lcg rng;
    rng.s = 2024;
    double worst = 0;
    int slowed = 0;
    for (int trial = 0; trial < 400; ++trial) {
        Body b = Body::solid_box({rng.range(0.03f, 0.4f), rng.range(0.03f, 0.4f), rng.range(0.03f, 0.4f)}, 1, {0, 0, 0}, rng.rotation());
        b.w = rng.vec(1).normalized() * rng.range(50, 900);  // up to 7.5 rad per step
        const double start = double(b.kinetic_energy());
        double peak = 0;
        for (int i = 0; i < 60; ++i) {
            b.integrate(1.0f / 120.0f, {}, Gyroscopic::Implicit);
            peak = std::max(peak, double(b.kinetic_energy()));
        }
        worst = std::max(worst, peak / start);
        slowed += double(b.kinetic_energy()) < 0.2 * start;
        CHECK(std::isfinite(b.w.x) && std::isfinite(b.w.y) && std::isfinite(b.w.z));
    }
    CHECK(worst <= 1.001);
    CHECK(slowed < 400);  // it is a guard, not a brake: not every spin is killed
}
