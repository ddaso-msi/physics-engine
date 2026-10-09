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
