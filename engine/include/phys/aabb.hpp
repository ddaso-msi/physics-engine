#pragma once
// Stage 7: axis-aligned bounding boxes, the cheap proxy the broad phase uses instead of real shapes.
#include "body.hpp"

#include <algorithm>

namespace phys {

struct AABB {
    Vec2 lo, hi;
};

// Touching boxes count as overlapping: better a spurious candidate than a missed contact.
inline bool overlap(const AABB& a, const AABB& b) {
    return a.lo.x <= b.hi.x && b.lo.x <= a.hi.x && a.lo.y <= b.hi.y && b.lo.y <= a.hi.y;
}

inline AABB merge(const AABB& a, const AABB& b) {
    return {{std::min(a.lo.x, b.lo.x), std::min(a.lo.y, b.lo.y)},
            {std::max(a.hi.x, b.hi.x), std::max(a.hi.y, b.hi.y)}};
}

// True when `inner` lies entirely inside `outer`.
inline bool contains(const AABB& outer, const AABB& inner) {
    return outer.lo.x <= inner.lo.x && outer.lo.y <= inner.lo.y && inner.hi.x <= outer.hi.x &&
           inner.hi.y <= outer.hi.y;
}

// Grow by `margin` on every side.
inline AABB fatten(const AABB& a, Real margin) {
    return {{a.lo.x - margin, a.lo.y - margin}, {a.hi.x + margin, a.hi.y + margin}};
}

// Half the true perimeter; the tree uses it as its "cost" of a box, since the chance a random query
// hits a box is proportional to its perimeter.
inline Real perimeter(const AABB& a) { return (a.hi.x - a.lo.x) + (a.hi.y - a.lo.y); }

// Tight world-space box around a body's shape at its current position and angle.
inline AABB compute_aabb(const Body& b) {
    if (b.shape.type == Shape::Type::Circle) {
        const Real r = b.shape.circle.radius;
        return {{b.pos.x - r, b.pos.y - r}, {b.pos.x + r, b.pos.y + r}};
    }
    const Polygon& p = b.shape.polygon;
    const Transform t = b.transform();
    Vec2 lo = apply(t, p.vertices[0]), hi = lo;
    for (int i = 1; i < p.count; ++i) {
        const Vec2 v = apply(t, p.vertices[i]);
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y)};
    }
    return {lo, hi};
}

}  // namespace phys
