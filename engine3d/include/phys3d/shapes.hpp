#pragma once
// Stage 10, step 3: collision shapes. A shape lives in its body's frame, centred on the centre of mass.
#include "math.hpp"

namespace phys3d {

struct Shape {
    enum class Type { Sphere, Box, Capsule };

    Type type = Type::Sphere;
    Real radius = 0;       // Sphere, Capsule
    Vec3 half_extents;     // Box: half the width along each of the body's own axes
    Real half_length = 0;  // Capsule: its axis runs along the body's own y, from -half_length to +half_length

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
    static Shape box(Vec3 half) {
        Shape s;
        s.type = Type::Box;
        s.half_extents = half;
        return s;
    }
};

}  // namespace phys3d
