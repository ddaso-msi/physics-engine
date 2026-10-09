#pragma once
// Stage 10, step 3: collision shapes. A shape lives in its body's frame, centred on the centre of mass.
#include "math.hpp"

namespace phys3d {

struct Shape {
    enum class Type { Sphere, Box };

    Type type = Type::Sphere;
    Real radius = 0;    // Sphere
    Vec3 half_extents;  // Box: half the width along each of the body's own axes

    static Shape sphere(Real r) {
        Shape s;
        s.type = Type::Sphere;
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
