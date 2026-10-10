#include <phys3d/body.hpp>

#include <cmath>

namespace phys3d {

Body Body::dynamic(Real m, const Mat3& inertia, Vec3 position, Quat orientation) {
    Body b;
    b.pos = position;
    b.q = orientation.normalized();
    b.mass = m;
    b.inv_mass = m > 0 ? 1 / m : 0;
    b.inertia_body = inertia;
    b.inv_inertia_body = inertia.inverse();  // zero matrix if the tensor is singular
    return b;
}

Body Body::solid_sphere(Real radius, Real density, Vec3 position) {
    const Real m = density * (static_cast<Real>(4) / 3) * kPi * radius * radius * radius;
    Body b = dynamic(m, sphere_inertia(m, radius), position);
    b.shape = Shape::sphere(radius);
    return b;
}

Body Body::solid_box(Vec3 half, Real density, Vec3 position, Quat orientation) {
    const Real m = density * 8 * half.x * half.y * half.z;
    Body b = dynamic(m, box_inertia(m, half), position, orientation);
    b.shape = Shape::box(half);
    return b;
}

Body Body::fixed(Vec3 position, Quat orientation) {
    Body b;
    b.type = BodyType::Static;
    b.pos = position;
    b.q = orientation.normalized();
    return b;
}

Body Body::fixed_box(Vec3 half, Vec3 position, Quat orientation) {
    Body b = fixed(position, orientation);
    b.shape = Shape::box(half);
    return b;
}

Body Body::fixed_sphere(Real radius, Vec3 position) {
    Body b = fixed(position);
    b.shape = Shape::sphere(radius);
    return b;
}

void Body::integrate_velocity(Real dt, Vec3 gravity, Gyroscopic mode) {
    if (type == BodyType::Dynamic) {
        vel += (force * inv_mass + gravity) * dt;

        const Mat3 rot = to_mat3(q);
        const Mat3 inv_i_world = to_world_frame(rot, inv_inertia_body);

        if (mode == Gyroscopic::Explicit) {
            // dw = I^-1 (torque - w x (I w)) dt, everything taken at the start of the step.
            const Vec3 l = to_world_frame(rot, inertia_body) * w;
            w += inv_i_world * (torque - cross(w, l)) * dt;
        } else {
            w += inv_i_world * torque * dt;
        }

        if (mode == Gyroscopic::Implicit) {
            // Work in the body frame, where the inertia tensor I is constant. With no torque,
            //     I (w1 - w0) + dt * w1 x (I w1) = 0
            // must hold for the END-of-step angular velocity w1. That is nonlinear in w1, so take one
            // Newton step from w0:  f(w) = I (w - w0) + dt w x (I w),  J = df/dw = I + dt ([w]x I - [I w]x),
            // w1 = w0 - J^-1 f(w0).  At w0 the first term of f vanishes, leaving f(w0) = dt w0 x (I w0).
            const Vec3 wb = rot.transposed() * w;
            const Vec3 iw = inertia_body * wb;
            const Vec3 f = cross(wb, iw) * dt;
            const Mat3 jac = inertia_body + (Mat3::skew(wb) * inertia_body - Mat3::skew(iw)) * dt;
            Vec3 wb1 = wb - jac.inverse() * f;

            // One Newton iteration is accurate while the body turns a modest angle per step. At extreme
            // spins (several radians per step, as after a hard corner impact) it can overshoot badly and
            // hand back far more spin than it was given. The gyroscopic term is a torque perpendicular to w:
            // it does no work, so the true answer has exactly the rotational energy we started with. Never
            // let the step exceed that; scale the result back onto it if it does.
            const Real energy0 = dot(wb, iw), energy1 = dot(wb1, inertia_body * wb1);
            if (!(energy1 <= energy0)) wb1 = energy1 > 0 && std::isfinite(energy1) ? wb1 * std::sqrt(energy0 / energy1) : wb;
            w = rot * wb1;
        }
    }
    force = {};
    torque = {};
}

void Body::integrate_position(Real dt) {
    if (type == BodyType::Static) return;
    pos += vel * dt;
    q = integrate_orientation(q, w, dt);
}

}  // namespace phys3d
