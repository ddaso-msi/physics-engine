#include "test.hpp"
#include <phys/integrate.hpp>
#include <phys/timestep.hpp>

using namespace phys;

namespace {

Vec2 gravity(const Particle&) { return {0, -9.81f}; }

// Unit-frequency harmonic oscillator, x'' = -x. Energy E = (v^2 + x^2)/2, exactly 0.5 at x=1, v=0.
Vec2 spring(const Particle& p) { return -p.pos; }
double energy(const Particle& p) { return 0.5 * (double(p.vel.length_sq()) + double(p.pos.length_sq())); }

Particle run_oscillator(Integrator kind, Real dt, int steps, double* max_energy = nullptr) {
    Particle p{{1, 0}, {0, 0}};
    for (int i = 0; i < steps; ++i) {
        step(p, dt, kind, spring);
        if (max_energy && energy(p) > *max_energy) *max_energy = energy(p);
    }
    return p;
}

}  // namespace

// With constant acceleration each scheme has a closed form after n steps of size dt.
// Exact answer is y = g*t^2/2 with t = n*dt.
TEST(free_fall_closed_forms) {
    const int n = 60;
    const Real dt = 1.0f / 60.0f;
    const double g = -9.81, t = n * double(dt);

    Particle e{{0, 0}, {0, 0}}, s{{0, 0}, {0, 0}}, v{{0, 0}, {0, 0}};
    for (int i = 0; i < n; ++i) {
        step(e, dt, Integrator::ExplicitEuler, gravity);
        step(s, dt, Integrator::SemiImplicitEuler, gravity);
        step(v, dt, Integrator::VelocityVerlet, gravity);
    }
    // explicit Euler lags the exact answer:       y = g*dt^2 * n(n-1)/2
    CHECK_NEAR(e.pos.y, g * double(dt) * dt * n * (n - 1) / 2, 1e-3);
    // semi-implicit Euler leads it:               y = g*dt^2 * n(n+1)/2
    CHECK_NEAR(s.pos.y, g * double(dt) * dt * n * (n + 1) / 2, 1e-3);
    // Verlet is exact for constant acceleration:  y = g*t^2/2
    CHECK_NEAR(v.pos.y, 0.5 * g * t * t, 1e-3);
    CHECK_NEAR(v.vel.y, g * t, 1e-3);
    // and the two Eulers straddle the truth
    CHECK(e.pos.y > 0.5 * g * t * t);   // less negative: fell too little
    CHECK(s.pos.y < 0.5 * g * t * t);   // fell too far
}

// 100 periods of the oscillator. Explicit Euler multiplies energy by (1 + dt^2) every step, so it
// grows like e^(2*pi*100*dt) ~ 535x here. The other two stay bounded near the true value 0.5.
TEST(oscillator_energy_drift) {
    const Real dt = 0.01f;
    const int steps = static_cast<int>(100 * 2 * 3.14159265 / double(dt));

    double e_explicit = 0, e_semi = 0, e_verlet = 0;
    double max_semi = 0, max_verlet = 0;
    e_explicit = energy(run_oscillator(Integrator::ExplicitEuler, dt, steps));
    e_semi = energy(run_oscillator(Integrator::SemiImplicitEuler, dt, steps, &max_semi));
    e_verlet = energy(run_oscillator(Integrator::VelocityVerlet, dt, steps, &max_verlet));

    CHECK(e_explicit > 100 * 0.5);
    CHECK_NEAR(e_semi, 0.5, 0.01);
    CHECK_NEAR(max_semi, 0.5, 0.01);       // bounded for the whole run, not just at the end
    CHECK_NEAR(e_verlet, 0.5, 1e-3);
    CHECK_NEAR(max_verlet, 0.5, 1e-3);
}

// Halving dt should cut the position error by ~2x for semi-implicit Euler (first order) and ~4x for
// Verlet (second order), measured against the exact solution x(t) = cos(t) after one second.
TEST(convergence_order) {
    auto error = [](Integrator kind, Real dt) {
        int steps = static_cast<int>(1.0 / double(dt) + 0.5);
        Particle p = run_oscillator(kind, dt, steps);
        return std::fabs(double(p.pos.x) - std::cos(1.0));
    };
    double semi_ratio = error(Integrator::SemiImplicitEuler, 0.02f) / error(Integrator::SemiImplicitEuler, 0.01f);
    double verlet_ratio = error(Integrator::VelocityVerlet, 0.02f) / error(Integrator::VelocityVerlet, 0.01f);
    CHECK_NEAR(semi_ratio, 2.0, 0.3);
    CHECK_NEAR(verlet_ratio, 4.0, 0.6);
}

// dt = 1/64 and these frame times are exact in binary floating point, so the counts are exact.
TEST(fixed_timestep_accumulator) {
    FixedTimestep ts(1.0f / 64.0f);
    int calls = 0;
    auto count = [&](Real dt) { CHECK_NEAR(dt, 1.0 / 64.0, 1e-9); ++calls; };

    CHECK(ts.advance(1.0f / 32.0f, count) == 2);        // two whole steps, nothing left over
    CHECK_NEAR(ts.alpha(), 0, 1e-9);

    CHECK(ts.advance(1.5f / 64.0f, count) == 1);        // 1.5 steps: one now, half carried over
    CHECK_NEAR(ts.alpha(), 0.5, 1e-9);
    CHECK(ts.advance(1.5f / 64.0f, count) == 2);        // carry + 1.5 = 2 whole steps
    CHECK_NEAR(ts.alpha(), 0, 1e-9);
    CHECK(calls == 5);

    CHECK(ts.advance(0.001f / 64.0f, count) == 0);      // too small to step: just banked
    CHECK(calls == 5);
}

TEST(fixed_timestep_caps_catch_up) {
    FixedTimestep ts(1.0f / 64.0f, 0.25f);
    int calls = 0;
    int n = ts.advance(10.0f, [&](Real) { ++calls; });  // a 10 s stall must not demand 640 steps
    CHECK(n == 16);                                     // 0.25 s / (1/64 s)
    CHECK(calls == 16);
}
