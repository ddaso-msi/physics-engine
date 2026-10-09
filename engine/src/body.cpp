#include <phys/body.hpp>

namespace phys {

Body::Body(const Shape& s, Vec2 position, Real a, BodyType t, Real density)
    : type(t), shape(s), pos(position) {
    set_angle(a);
    if (type == BodyType::Dynamic) {
        MassData md = compute_mass(shape, density);
        mass = md.mass;
        inertia = md.inertia;
        inv_mass = mass > 0 ? 1 / mass : 0;
        inv_inertia = inertia > 0 ? 1 / inertia : 0;
    }
}

Transform Body::transform() const {
    Transform t;
    t.p = pos;
    t.q = q;
    return t;
}

void Body::set_angle(Real a) {
    // Keep the stored angle in (-pi, pi] so float precision does not degrade as a body keeps spinning.
    angle = a - 2 * kPi * std::round(a / (2 * kPi));
    q = Rot(angle);
}

Real Body::inscribed_radius() const {
    if (shape.type == Shape::Type::Circle) return shape.circle.radius;
    // The centre of mass is the polygon's origin, so the distance to a face is dot(normal, any vertex on it).
    const Polygon& p = shape.polygon;
    Real r = dot(p.normals[0], p.vertices[0]);
    for (int i = 1; i < p.count; ++i) r = std::min(r, dot(p.normals[i], p.vertices[i]));
    return r;
}

void Body::apply_force_at(Vec2 f, Vec2 world_point) {
    force += f;
    torque += cross(world_point - pos, f);
}

void Body::apply_impulse_at(Vec2 j, Vec2 world_point) {
    vel += j * inv_mass;
    w += inv_inertia * cross(world_point - pos, j);
}

Real Body::kinetic_energy() const {
    return static_cast<Real>(0.5) * (mass * vel.length_sq() + inertia * w * w);
}

bool Body::contains(Vec2 world_point) const {
    const Vec2 local = inv_rotate(q, world_point - pos);
    if (shape.type == Shape::Type::Circle) return local.length_sq() <= shape.circle.radius * shape.circle.radius;
    const Polygon& poly = shape.polygon;
    for (int i = 0; i < poly.count; ++i)
        if (dot(poly.normals[i], local - poly.vertices[i]) > 0) return false;
    return true;
}

void Body::integrate_velocity(Real dt, Vec2 gravity) {
    if (type == BodyType::Dynamic) {
        vel += (force * inv_mass + gravity) * dt;
        w += torque * inv_inertia * dt;
    }
    force = {};
    torque = 0;
}

void Body::integrate_position(Real dt) {
    if (type == BodyType::Static) return;
    pos += vel * dt;
    set_angle(angle + w * dt);
}

void Body::integrate(Real dt, Vec2 gravity) {
    integrate_velocity(dt, gravity);
    integrate_position(dt);
}

}  // namespace phys
