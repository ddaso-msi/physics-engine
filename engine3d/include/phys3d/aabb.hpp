#pragma once
// Stage 10, step 6: axis-aligned bounding boxes in 3D, the cheap stand-in the broad phase compares.
#include "body.hpp"

#include <algorithm>

namespace phys3d {

struct AABB {
    Vec3 lo, hi;
};

// Touching boxes count as overlapping: better a spurious candidate than a missed contact.
inline bool overlap(const AABB& a, const AABB& b) {
    return a.lo.x <= b.hi.x && b.lo.x <= a.hi.x && a.lo.y <= b.hi.y && b.lo.y <= a.hi.y && a.lo.z <= b.hi.z && b.lo.z <= a.hi.z;
}

inline AABB merge(const AABB& a, const AABB& b) {
    return {{std::min(a.lo.x, b.lo.x), std::min(a.lo.y, b.lo.y), std::min(a.lo.z, b.lo.z)},
            {std::max(a.hi.x, b.hi.x), std::max(a.hi.y, b.hi.y), std::max(a.hi.z, b.hi.z)}};
}

// True when `inner` lies entirely inside `outer`.
inline bool contains(const AABB& outer, const AABB& inner) {
    return outer.lo.x <= inner.lo.x && outer.lo.y <= inner.lo.y && outer.lo.z <= inner.lo.z && inner.hi.x <= outer.hi.x &&
           inner.hi.y <= outer.hi.y && inner.hi.z <= outer.hi.z;
}

// Grow by `margin` on every side.
inline AABB fatten(const AABB& a, Real margin) {
    return {a.lo - Vec3{margin, margin, margin}, a.hi + Vec3{margin, margin, margin}};
}

// Half the surface area. The tree uses it as the "cost" of a box: the chance that a random query box or
// ray touches a box is proportional to its surface area (the 3D counterpart of the perimeter used in 2D).
inline Real surface_area(const AABB& a) {
    const Vec3 d = a.hi - a.lo;
    return d.x * d.y + d.y * d.z + d.z * d.x;
}

// Tight world-space box around a body's shape. For a turned box, the extent along each world axis is the
// sum of the box's three half-extents, each scaled by how much its axis leans along that world axis:
// extent_i = sum_j |R_ij| h_j.
inline AABB compute_aabb(const Body& b) {
    if (b.shape.type != Shape::Type::Box) {
        // A sphere, or a capsule: the box around its axis (a single point for a sphere), grown by the radius.
        const Vec3 r{b.shape.radius, b.shape.radius, b.shape.radius};
        const Vec3 h = rotate(b.q, {0, b.shape.half_length, 0});
        const Vec3 reach{std::fabs(h.x), std::fabs(h.y), std::fabs(h.z)};
        return {b.pos - reach - r, b.pos + reach + r};
    }
    const Mat3 rot = to_mat3(b.q);
    const Vec3 h = b.shape.half_extents;
    Vec3 e;
    for (int i = 0; i < 3; ++i) e[i] = std::fabs(rot.at(i, 0)) * h.x + std::fabs(rot.at(i, 1)) * h.y + std::fabs(rot.at(i, 2)) * h.z;
    return {b.pos - e, b.pos + e};
}

}  // namespace phys3d
