#pragma once
// Stage 3: rigid body state and its integration.
//
// pos is the centre of mass. Forces and torques accumulate until integrate() consumes them.
// Gravity is applied as an acceleration, so it does not depend on mass.
#include "math.hpp"
#include "shapes.hpp"

#include <cstdint>

namespace phys {

enum class BodyType { Static, Dynamic };

struct Body {
    BodyType type = BodyType::Dynamic;
    Shape shape;

    Vec2 pos;
    Real angle = 0;
    Rot q;

    Vec2 vel;
    Real w = 0;  // angular velocity, rad/s, counter-clockwise positive

    Vec2 force;
    Real torque = 0;

    // Inverse mass/inertia of 0 means "infinite": the body cannot be accelerated or spun.
    Real mass = 0, inv_mass = 0;
    Real inertia = 0, inv_inertia = 0;

    // Surface material. A contact between two bodies mixes them (see solver.cpp).
    Real restitution = static_cast<Real>(0.2);  // 0 = no bounce, 1 = perfectly elastic
    Real friction = static_cast<Real>(0.5);     // Coulomb coefficient

    // Collision filtering. `category` says which group(s) this body belongs to, `mask` which groups
    // it is willing to collide with. Two bodies collide only if EACH one's mask admits the other's
    // category, so a one-sided mask cannot make a body push another that cannot push back.
    std::uint16_t category = 0x0001;
    std::uint16_t mask = 0xFFFF;

    Body() = default;
    Body(const Shape& shape, Vec2 position, Real angle, BodyType type = BodyType::Dynamic,
         Real density = 1);

    Transform transform() const;
    void set_angle(Real a);

    void apply_force(Vec2 f) { force += f; }
    void apply_torque(Real t) { torque += t; }
    // A force away from the centre of mass also produces torque r x F.
    void apply_force_at(Vec2 f, Vec2 world_point);
    // Instant change of momentum: dv = j/m, dw = (r x j)/I.
    void apply_impulse_at(Vec2 j, Vec2 world_point);

    Vec2 velocity_at(Vec2 world_point) const { return vel + cross(w, world_point - pos); }
    Real kinetic_energy() const;
    bool contains(Vec2 world_point) const;

    // Semi-implicit Euler (see integrate.hpp), split in two so a contact solver can run between the
    // halves: velocities first (consumes and clears forces), then positions from those velocities.
    void integrate_velocity(Real dt, Vec2 gravity);
    void integrate_position(Real dt);
    // Both halves back to back, for bodies that are not in a contact-solving world.
    void integrate(Real dt, Vec2 gravity);
};

inline bool should_collide(const Body& a, const Body& b) {
    return (a.mask & b.category) != 0 && (b.mask & a.category) != 0;
}

}  // namespace phys
