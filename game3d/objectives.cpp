#include "objectives.hpp"

#include <phys3d/aabb.hpp>

#include <algorithm>

namespace puzzle3d {

int count_below(const World& w, std::span<const int> bodies, Real y) {
    int n = 0;
    for (int i : bodies)
        if (w.bodies[static_cast<size_t>(i)].pos.y < y) ++n;
    return n;
}

bool inside(const World& w, int body, Vec3 lo, Vec3 hi) {
    const Vec3 p = w.bodies[static_cast<size_t>(body)].pos;
    return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z;
}

bool all_slow(const World& w, std::span<const int> bodies, Real linear, Real angular) {
    for (int i : bodies) {
        const Body& b = w.bodies[static_cast<size_t>(i)];
        if (b.vel.length() > linear || b.w.length() > angular) return false;
    }
    return true;
}

Real top_of(const World& w, std::span<const int> bodies) {
    Real top = kGroundY;
    for (int i : bodies) top = std::max(top, phys3d::compute_aabb(w.bodies[static_cast<size_t>(i)]).hi.y);
    return top;
}

}  // namespace puzzle3d
