#pragma once
// Stage 2: numerical integration of a point mass.
//
// We solve  x'' = a(x, v)  by stepping in fixed increments dt. The three schemes differ only in the
// ORDER of the position/velocity updates, and that order decides whether energy blows up, stays
// bounded, or is conserved to second order.
#include "math.hpp"

namespace phys {

struct Particle {
    Vec2 pos;
    Vec2 vel;
};

enum class Integrator {
    // x += v*dt ; v += a*dt          (uses the OLD velocity for position). Gains energy: unstable.
    ExplicitEuler,
    // v += a*dt ; x += v*dt          (uses the NEW velocity for position). Symplectic: energy stays
    // bounded. What most game physics engines (Box2D included) use.
    SemiImplicitEuler,
    // x += v*dt + a*dt^2/2 ; v += (a_old + a_new)*dt/2.  Second-order accurate, symplectic, exact for
    // constant acceleration. Needs the acceleration evaluated twice per step.
    VelocityVerlet,
};

// `accel` is any callable `Vec2(const Particle&)` returning the acceleration for that state.
template <class Accel>
void step(Particle& p, Real dt, Integrator kind, Accel&& accel) {
    switch (kind) {
        case Integrator::ExplicitEuler: {
            Vec2 a = accel(p);
            p.pos += p.vel * dt;
            p.vel += a * dt;
            break;
        }
        case Integrator::SemiImplicitEuler: {
            Vec2 a = accel(p);
            p.vel += a * dt;
            p.pos += p.vel * dt;
            break;
        }
        case Integrator::VelocityVerlet: {
            Vec2 a0 = accel(p);
            p.pos += p.vel * dt + a0 * (dt * dt * static_cast<Real>(0.5));
            Vec2 a1 = accel(p);  // at the new position (velocity still old: fine for position-only forces)
            p.vel += (a0 + a1) * (dt * static_cast<Real>(0.5));
            break;
        }
    }
}

}  // namespace phys
