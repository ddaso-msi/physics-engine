#pragma once
// A vector in double precision, for the code that has to work below what single precision can resolve: the
// hull builder (step 12) and GJK and EPA (step 18). Internal to the engine; the engine's Real is used for
// what goes in and what comes out.
#include <phys3d/math.hpp>

#include <cmath>

namespace phys3d {

struct Vec3d {
    double x = 0, y = 0, z = 0;
    Vec3d() = default;
    Vec3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    explicit Vec3d(Vec3 v) : x(v.x), y(v.y), z(v.z) {}

    double length_sq() const { return x * x + y * y + z * z; }
    double length() const { return std::sqrt(length_sq()); }
    Vec3d normalized() const {
        const double len = length();
        return len > 0 ? Vec3d{x / len, y / len, z / len} : Vec3d{};
    }
    Vec3d& operator+=(Vec3d b) { return x += b.x, y += b.y, z += b.z, *this; }
};

inline Vec3d operator+(Vec3d a, Vec3d b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3d operator-(Vec3d a, Vec3d b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3d operator-(Vec3d a) { return {-a.x, -a.y, -a.z}; }
inline Vec3d operator*(Vec3d a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3d operator/(Vec3d a, double s) { return {a.x / s, a.y / s, a.z / s}; }
inline double dot(Vec3d a, Vec3d b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3d cross(Vec3d a, Vec3d b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline double distance(Vec3d a, Vec3d b) { return (a - b).length(); }
inline Vec3 to_real(Vec3d v) { return {static_cast<Real>(v.x), static_cast<Real>(v.y), static_cast<Real>(v.z)}; }

// The same formula as rotate() for Vec3, so for a quaternion that is not quite unit length the two still
// agree (to rounding); and rotate_d(conjugate(q), .) is exactly the transpose of rotate_d(q, .).
inline Vec3d rotate_d(const Quat& q, Vec3d v) {
    const Vec3d u{q.x, q.y, q.z};
    const Vec3d t = cross(u, v) * 2.0;
    return v + t * static_cast<double>(q.w) + cross(u, t);
}

}  // namespace phys3d
