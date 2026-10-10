#pragma once
// Stage 10, step 3: collision shapes. A shape lives in its body's frame, centred on the centre of mass.
#include "hull.hpp"
#include "math.hpp"

#include <algorithm>
#include <memory>

namespace phys3d {

struct Shape {
    enum class Type { Sphere, Box, Capsule, Hull };

    Type type = Type::Sphere;
    Real radius = 0;       // Sphere, Capsule
    Vec3 half_extents;     // Box: half the width along each of the body's own axes
    Real half_length = 0;  // Capsule: its axis runs along the body's own y, from -half_length to +half_length
    // Hull: a convex polyhedron, its corners measured from its centre of mass. Shared, not copied, when a
    // body is: a hundred identical crates point at one hull.
    std::shared_ptr<const Hull> hull;
    Real hull_inscribed = 0;  // the largest sphere about the centre that fits inside: the distance to the nearest face
    Real hull_bounding = 0;   // the smallest that contains it: the distance to the farthest corner

    static Shape sphere(Real r) {
        Shape s;
        s.type = Type::Sphere;
        s.radius = r;
        return s;
    }
    // Every point within `r` of a segment: a cylinder with a hemisphere on each end. Its overall length is
    // 2 * (half_len + r).
    static Shape capsule(Real half_len, Real r) {
        Shape s;
        s.type = Type::Capsule;
        s.half_length = half_len;
        s.radius = r;
        return s;
    }
    static Shape convex(std::shared_ptr<const Hull> centred) {
        Shape s;
        s.type = Type::Hull;
        s.hull_inscribed = centred->faces.empty() ? Real(0) : centred->faces[0].offset;
        for (const Hull::Face& f : centred->faces) s.hull_inscribed = std::min(s.hull_inscribed, f.offset);
        for (const Vec3& v : centred->vertices) s.hull_bounding = std::max(s.hull_bounding, v.length());
        s.hull = std::move(centred);
        return s;
    }
    static Shape box(Vec3 half) {
        Shape s;
        s.type = Type::Box;
        s.half_extents = half;
        return s;
    }
};

}  // namespace phys3d
