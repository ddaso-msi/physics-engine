#pragma once
// Stage 10, step 1: 3D inertia tensors.
//
// In 2D a body's resistance to spinning is one number. In 3D it is a symmetric 3x3 matrix I that turns an
// angular velocity into an angular momentum, L = I w, and in general L is NOT parallel to w. The tensor is
// constant in the body's own frame, but the world-frame tensor changes as the body turns:  I_world = R I R^T.
// Solvers need the INVERSE (acceleration = I^-1 * torque), so the engine stores inverse_inertia in the body
// frame and rotates that.
#include "math.hpp"

namespace phys3d {

// Uniform solid sphere: I = 2/5 m r^2 about any axis through the centre.
inline Mat3 sphere_inertia(Real mass, Real radius) {
    const Real i = static_cast<Real>(0.4) * mass * radius * radius;
    return Mat3::diagonal(i, i, i);
}

// Uniform solid box with HALF extents (a, b, c), about its centre and aligned with its axes:
// Ixx = m (b^2 + c^2) / 3 and so on (the usual m (w^2 + h^2) / 12 with full widths w = 2b, h = 2c).
inline Mat3 box_inertia(Real mass, Vec3 half) {
    const Real a2 = half.x * half.x, b2 = half.y * half.y, c2 = half.z * half.z;
    return Mat3::diagonal(mass * (b2 + c2) / 3, mass * (a2 + c2) / 3, mass * (a2 + b2) / 3);
}

// Uniform solid capsule with its axis along y: a cylinder of length 2 * half_length and radius r, plus a
// hemisphere on each end. Split the mass by volume into the cylinder's share mc and the two hemispheres'
// (together one sphere's) ms. Then
//   about the axis:   Iyy = mc r^2 / 2 + ms 2 r^2 / 5                    (a cylinder and a sphere)
//   across the axis:  Ixx = Izz = mc (L^2 / 12 + r^2 / 4) + ms (2 r^2 / 5 + h^2 + 3 h r / 4),   L = 2 h.
// The last bracket is the two hemispheres moved out to the ends: a hemisphere's centre of mass is 3r/8
// from its flat face and its own transverse inertia about that point is 83/320 m r^2, so by the parallel
// axis theorem each half contributes 83/320 r^2 + (h + 3r/8)^2 per unit mass, which multiplies out to the above.
inline Mat3 capsule_inertia(Real mass, Real half_length, Real radius) {
    const Real cylinder_volume = 2 * half_length * radius * radius, sphere_volume = (static_cast<Real>(4) / 3) * radius * radius * radius;  // both over pi
    const Real mc = mass * cylinder_volume / (cylinder_volume + sphere_volume), ms = mass - mc;
    const Real r2 = radius * radius, h = half_length;
    const Real about_axis = mc * r2 / 2 + ms * static_cast<Real>(0.4) * r2;
    const Real across = mc * (h * h / 3 + r2 / 4) + ms * (static_cast<Real>(0.4) * r2 + h * h + static_cast<Real>(0.75) * h * radius);
    return Mat3::diagonal(across, about_axis, across);
}
inline Real capsule_volume(Real half_length, Real radius) {
    return kPi * radius * radius * (2 * half_length + (static_cast<Real>(4) / 3) * radius);
}

// Inertia about a point offset by `d` from the centre of mass (parallel axis theorem):
// I' = I + m (|d|^2 E - d d^T), with E the identity.
inline Mat3 parallel_axis(const Mat3& inertia_at_centre, Real mass, Vec3 d) {
    return inertia_at_centre + (Mat3::identity() * d.length_sq() - Mat3::outer(d, d)) * mass;
}

// Rotate a body-frame tensor into the world frame for a body whose orientation is `rotation`. Works for
// the inertia tensor and, equally, for its inverse.
inline Mat3 to_world_frame(const Mat3& rotation, const Mat3& body_frame_tensor) {
    return rotation * body_frame_tensor * rotation.transposed();
}

}  // namespace phys3d
