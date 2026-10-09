#pragma once
// Stage 3: collision shapes and the mass properties derived from them.
//
// Convention: a shape lives in BODY space with its centre of mass at the origin. Polygons are
// recentred on their centroid when built, so a body's position is always its centre of mass and
// the inertia we compute is already about that point.
#include "math.hpp"

#include <optional>
#include <span>

namespace phys {

struct MassData {
    Real mass = 0;
    Real inertia = 0;  // about the centre of mass (the body origin)
};

struct Circle {
    Real radius = 0;
};

// Convex polygon, counter-clockwise, at most kMaxVertices. normals[i] is the outward unit normal of
// the edge from vertices[i] to vertices[i+1].
struct Polygon {
    static constexpr int kMaxVertices = 8;

    int count = 0;
    Vec2 vertices[kMaxVertices];
    Vec2 normals[kMaxVertices];
    // Where the centroid was in the caller's coordinates, i.e. how far the points were shifted.
    Vec2 centroid_shift;

    // Accepts either winding. Fails (nullopt) for <3 or >kMaxVertices points, zero area, or any
    // non-convex / collinear corner.
    static std::optional<Polygon> from_points(std::span<const Vec2> points);
    // Axis-aligned box centred on the origin with the given half extents.
    static Polygon box(Real half_width, Real half_height);
};

struct Shape {
    enum class Type { Circle, Polygon };

    Type type = Type::Circle;
    Circle circle;
    Polygon polygon;

    static Shape make_circle(Real radius);
    static Shape make_polygon(const Polygon& polygon);
};

// Mass and inertia of a uniform-density shape (2D, so density is mass per unit area).
//   circle:  m = rho*pi*r^2,  I = m*r^2/2
//   polygon: split into triangles (origin, v_i, v_i+1); each contributes
//            area = D/2,  I = rho * D/12 * (e1.e1 + e1.e2 + e2.e2)   with D = cross(e1, e2).
MassData compute_mass(const Shape& shape, Real density);

}  // namespace phys
