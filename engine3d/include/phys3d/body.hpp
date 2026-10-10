#pragma once
// Stage 10, step 2: 3D rigid body state and its integration.
//
// Like the 2D Body: `pos` is the centre of mass, forces and torques accumulate until integrate_velocity
// consumes them, gravity is an acceleration, and inverse mass / inverse inertia of zero mean "immovable".
// What is new in 3D:
//   - orientation is a quaternion `q`, and angular velocity `w` is a vector (world frame, rad/s);
//   - inertia is a tensor. It is stored in the BODY frame (where it is constant) and rotated into the world
//     frame when needed: I_world = R I_body R^T;
//   - angular momentum L = I_world w is what a free body conserves, not w. Because I_world turns with the
//     body, w must change even with no torque at all:  dw/dt = I^-1 (torque - w x (I w)).
//     That last term is the gyroscopic term. It is why a thrown book tumbles.
#include "inertia.hpp"
#include "shapes.hpp"

#include <algorithm>

namespace phys3d {

enum class BodyType { Static, Dynamic };

// How integrate_velocity treats the gyroscopic term.
enum class Gyroscopic {
    Off,       // ignore it: w stays constant without torque. Cheap, and wrong for anything not ball-like.
    Explicit,  // use the w at the start of the step. Simple, but it GAINS energy and blows up on fast spins.
    Implicit,  // solve for the w at the end of the step (one Newton iteration in the body frame). Stable:
               // it can only lose a little energy. The default.
};

struct Body {
    BodyType type = BodyType::Dynamic;
    Shape shape;  // what it collides as (step 3)

    Vec3 pos;
    Quat q;

    Vec3 vel;
    Vec3 w;  // angular velocity, world frame: the body turns about w at |w| rad/s

    Vec3 force;
    Vec3 torque;

    Real mass = 0, inv_mass = 0;
    Mat3 inertia_body = Mat3::diagonal(0, 0, 0);      // about the centre of mass, in the body's own axes
    Mat3 inv_inertia_body = Mat3::diagonal(0, 0, 0);  // zero = cannot be spun

    Real restitution = static_cast<Real>(0.2);
    Real friction = static_cast<Real>(0.5);

    // Sleeping (see World::step). A body whose whole island has been nearly still for long enough stops
    // being simulated until something disturbs it. Do not move a sleeping body by hand without wake().
    bool awake = true;
    bool allow_sleep = true;  // false: this body keeps its whole island awake
    Real sleep_time = 0;      // how long it has been nearly still

    Body() = default;

    // A dynamic body with the given mass and body-frame inertia tensor.
    static Body dynamic(Real mass, const Mat3& inertia_body, Vec3 pos, Quat q = {});
    // Uniform solid shapes: mass, inertia and collision shape all follow from the geometry.
    static Body solid_sphere(Real radius, Real density, Vec3 pos);
    static Body solid_box(Vec3 half_extents, Real density, Vec3 pos, Quat q = {});
    // The capsule's axis is the body's own y axis.
    static Body solid_capsule(Real half_length, Real radius, Real density, Vec3 pos, Quat q = {});
    // An immovable body (give it a shape, or use the helpers below).
    static Body fixed(Vec3 pos, Quat q = {});
    static Body fixed_box(Vec3 half_extents, Vec3 pos, Quat q = {});
    static Body fixed_sphere(Real radius, Vec3 pos);
    static Body fixed_capsule(Real half_length, Real radius, Vec3 pos, Quat q = {});

    Transform transform() const { return {pos, q}; }

    // Dynamic and not asleep: the bodies that actually get simulated.
    bool is_active() const { return type == BodyType::Dynamic && awake; }
    void wake() {
        awake = true;
        sleep_time = 0;
    }
    // Radius of the largest sphere around the centre of mass that fits inside the shape, and of the
    // smallest one that contains it. A body that moves less than the first in a step cannot have passed
    // straight through anything; everything it can touch lies within the second.
    Real inscribed_radius() const {
        return shape.type == Shape::Type::Box ? std::min({shape.half_extents.x, shape.half_extents.y, shape.half_extents.z}) : shape.radius;
    }
    Real bounding_radius() const {
        return shape.type == Shape::Type::Box ? shape.half_extents.length() : shape.radius + shape.half_length;  // (a sphere's half_length is 0)
    }
    Mat3 inertia_world() const { return to_world_frame(to_mat3(q), inertia_body); }
    Mat3 inv_inertia_world() const { return to_world_frame(to_mat3(q), inv_inertia_body); }

    Vec3 angular_momentum() const { return inertia_world() * w; }
    Real kinetic_energy() const {
        return static_cast<Real>(0.5) * (mass * vel.length_sq() + dot(w, angular_momentum()));
    }
    // Velocity of the material point of the body currently at `world_point`.
    Vec3 velocity_at(Vec3 world_point) const { return vel + cross(w, world_point - pos); }

    void apply_force(Vec3 f) { force += f; }
    void apply_torque(Vec3 t) { torque += t; }
    // A force away from the centre of mass also produces the torque r x F.
    void apply_force_at(Vec3 f, Vec3 world_point) {
        force += f;
        torque += cross(world_point - pos, f);
    }
    // Instant change of momentum: dv = j / m, dw = I^-1 (r x j).
    void apply_impulse_at(Vec3 j, Vec3 world_point) {
        vel += j * inv_mass;
        w += inv_inertia_world() * cross(world_point - pos, j);
    }

    // Semi-implicit Euler in two halves, as in 2D, so a contact solver can run between them. The velocity
    // half consumes and clears the force and torque.
    void integrate_velocity(Real dt, Vec3 gravity, Gyroscopic mode = Gyroscopic::Implicit);
    void integrate_position(Real dt);
    void integrate(Real dt, Vec3 gravity, Gyroscopic mode = Gyroscopic::Implicit) {
        integrate_velocity(dt, gravity, mode);
        integrate_position(dt);
    }
};

}  // namespace phys3d
