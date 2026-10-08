#pragma once
// Stage 1: 2D math. Everything is plain structs with value semantics.
// Conventions: right-handed, y up, angles in radians, counter-clockwise positive.
#include <cmath>

namespace phys {

using Real = float;

constexpr Real kPi = static_cast<Real>(3.14159265358979323846);

struct Vec2 {
    Real x = 0, y = 0;

    constexpr Vec2() = default;
    constexpr Vec2(Real x_, Real y_) : x(x_), y(y_) {}

    constexpr Vec2 operator-() const { return {-x, -y}; }
    constexpr Vec2& operator+=(Vec2 o) { x += o.x; y += o.y; return *this; }
    constexpr Vec2& operator-=(Vec2 o) { x -= o.x; y -= o.y; return *this; }
    constexpr Vec2& operator*=(Real s) { x *= s; y *= s; return *this; }

    Real length() const { return std::sqrt(x * x + y * y); }
    constexpr Real length_sq() const { return x * x + y * y; }

    // Returns the zero vector for (near-)zero input rather than NaN.
    Vec2 normalized() const {
        Real len = length();
        return len > static_cast<Real>(1e-12) ? Vec2{x / len, y / len} : Vec2{};
    }
};

constexpr Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
constexpr Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
constexpr Vec2 operator*(Vec2 v, Real s) { return {v.x * s, v.y * s}; }
constexpr Vec2 operator*(Real s, Vec2 v) { return {v.x * s, v.y * s}; }

constexpr Real dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }

// 2D "cross products". The z-component of the 3D cross product of two vectors in the plane.
constexpr Real cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }
// Angular velocity w (about +z) applied to offset r gives linear velocity w x r.
constexpr Vec2 cross(Real w, Vec2 r) { return {-w * r.y, w * r.x}; }
constexpr Vec2 cross(Vec2 r, Real w) { return {w * r.y, -w * r.x}; }

// Rotate 90 degrees counter-clockwise.
constexpr Vec2 perp(Vec2 v) { return {-v.y, v.x}; }

inline Real distance(Vec2 a, Vec2 b) { return (a - b).length(); }

// A rotation stored as sin/cos so it is built once and applied cheaply.
struct Rot {
    Real s = 0, c = 1;

    Rot() = default;
    explicit Rot(Real angle) : s(std::sin(angle)), c(std::cos(angle)) {}

    Real angle() const { return std::atan2(s, c); }
    Vec2 x_axis() const { return {c, s}; }
    Vec2 y_axis() const { return {-s, c}; }
};

inline Vec2 rotate(Rot q, Vec2 v) { return {q.c * v.x - q.s * v.y, q.s * v.x + q.c * v.y}; }
inline Vec2 inv_rotate(Rot q, Vec2 v) { return {q.c * v.x + q.s * v.y, -q.s * v.x + q.c * v.y}; }

// 2x2 matrix, column-major: columns ex and ey.
struct Mat2 {
    Vec2 ex{1, 0}, ey{0, 1};

    constexpr Mat2() = default;
    constexpr Mat2(Vec2 c1, Vec2 c2) : ex(c1), ey(c2) {}

    constexpr Real det() const { return ex.x * ey.y - ey.x * ex.y; }

    // Caller must ensure det() != 0; returns zero matrix when singular.
    constexpr Mat2 inverse() const {
        Real d = det();
        if (d != 0) d = 1 / d;
        return {{d * ey.y, -d * ex.y}, {-d * ey.x, d * ex.x}};
    }

    constexpr Mat2 transposed() const { return {{ex.x, ey.x}, {ex.y, ey.y}}; }
};

constexpr Vec2 operator*(const Mat2& m, Vec2 v) {
    return {m.ex.x * v.x + m.ey.x * v.y, m.ex.y * v.x + m.ey.y * v.y};
}

// Rigid transform: rotate then translate.
struct Transform {
    Vec2 p;
    Rot q;

    Transform() = default;
    Transform(Vec2 pos, Real angle) : p(pos), q(angle) {}
};

inline Vec2 apply(const Transform& t, Vec2 local) { return rotate(t.q, local) + t.p; }
inline Vec2 apply_inv(const Transform& t, Vec2 world) { return inv_rotate(t.q, world - t.p); }

}  // namespace phys
