#include <phys/shapes.hpp>

#include <algorithm>

namespace phys {

namespace {
constexpr Real kMinArea = static_cast<Real>(1e-6);
}

std::optional<Polygon> Polygon::from_points(std::span<const Vec2> points) {
    const int n = static_cast<int>(points.size());
    if (n < 3 || n > kMaxVertices) return std::nullopt;

    Vec2 pts[kMaxVertices];
    std::copy(points.begin(), points.end(), pts);

    // Signed area (positive = counter-clockwise) and centroid, as a triangle fan around the vertex
    // average. Using a reference point inside the polygon keeps the sums well-conditioned.
    Vec2 ref;
    for (int i = 0; i < n; ++i) ref += pts[i];
    ref *= 1 / static_cast<Real>(n);

    Real area = 0;
    Vec2 centroid;
    for (int i = 0; i < n; ++i) {
        Vec2 e1 = pts[i] - ref, e2 = pts[(i + 1) % n] - ref;
        Real tri = cross(e1, e2) * static_cast<Real>(0.5);
        area += tri;
        centroid += (e1 + e2) * (tri / 3);
    }
    if (std::fabs(area) < kMinArea) return std::nullopt;
    if (area < 0) std::reverse(pts, pts + n);  // normalise to counter-clockwise
    centroid *= 1 / std::fabs(area);
    centroid += ref;

    Polygon poly;
    poly.count = n;
    poly.centroid_shift = centroid;
    for (int i = 0; i < n; ++i) poly.vertices[i] = pts[i] - centroid;

    // Convex and strictly so: every corner must turn left.
    for (int i = 0; i < n; ++i) {
        Vec2 a = poly.vertices[(i + 1) % n] - poly.vertices[i];
        Vec2 b = poly.vertices[(i + 2) % n] - poly.vertices[(i + 1) % n];
        if (cross(a, b) <= 0) return std::nullopt;
    }
    for (int i = 0; i < n; ++i) {
        Vec2 edge = poly.vertices[(i + 1) % n] - poly.vertices[i];
        poly.normals[i] = Vec2{edge.y, -edge.x}.normalized();  // edge rotated 90 degrees clockwise
    }
    return poly;
}

Polygon Polygon::box(Real hx, Real hy) {
    const Vec2 pts[4] = {{-hx, -hy}, {hx, -hy}, {hx, hy}, {-hx, hy}};
    return *from_points(pts);
}

Shape Shape::make_circle(Real radius) {
    Shape s;
    s.type = Type::Circle;
    s.circle.radius = radius;
    return s;
}

Shape Shape::make_polygon(const Polygon& polygon) {
    Shape s;
    s.type = Type::Polygon;
    s.polygon = polygon;
    return s;
}

MassData compute_mass(const Shape& shape, Real density) {
    MassData md;
    if (shape.type == Shape::Type::Circle) {
        const Real r = shape.circle.radius;
        md.mass = density * kPi * r * r;
        md.inertia = static_cast<Real>(0.5) * md.mass * r * r;
        return md;
    }
    const Polygon& p = shape.polygon;
    Real area = 0, polar = 0;
    for (int i = 0; i < p.count; ++i) {
        Vec2 e1 = p.vertices[i], e2 = p.vertices[(i + 1) % p.count];
        Real d = cross(e1, e2);
        area += static_cast<Real>(0.5) * d;
        polar += d * (dot(e1, e1) + dot(e1, e2) + dot(e2, e2)) / 12;
    }
    md.mass = density * area;
    md.inertia = density * polar;
    return md;
}

}  // namespace phys
