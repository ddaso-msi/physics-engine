#pragma once
// Stage 10, step 1: 3D math. Vec3, Mat3, Quat and a rigid Transform.
//
// Conventions (the same ones the 2D engine uses, extended):
//   - right-handed axes, y up; a positive rotation about an axis turns counter-clockwise when you look
//     at the origin from the positive end of that axis (x -> y -> z -> x);
//   - Mat3 is stored by COLUMNS and multiplies column vectors:  (M v) = cx*v.x + cy*v.y + cz*v.z;
//   - Quat is (x, y, z, w) with w the scalar part. A unit quaternion for a rotation of `angle` about unit
//     `axis` is (axis*sin(angle/2), cos(angle/2)). q and -q are the same rotation.
// Real is the 2D engine's Real (float), so the two engines flip precision together.
#include <phys/math.hpp>

#include <cmath>

namespace phys3d {

using phys::kPi;
using phys::Real;

struct Vec3 {
    Real x = 0, y = 0, z = 0;

    constexpr Vec3() = default;
    constexpr Vec3(Real x_, Real y_, Real z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator-() const { return {-x, -y, -z}; }
    constexpr Vec3& operator+=(Vec3 o) { x += o.x; y += o.y; z += o.z; return *this; }
    constexpr Vec3& operator-=(Vec3 o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    constexpr Vec3& operator*=(Real s) { x *= s; y *= s; z *= s; return *this; }

    // Component by index 0, 1, 2 (for loops over axes).
    constexpr Real operator[](int i) const { return i == 0 ? x : i == 1 ? y : z; }
    constexpr Real& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }

    constexpr Real length_sq() const { return x * x + y * y + z * z; }
    Real length() const { return std::sqrt(length_sq()); }
    // The zero vector for (near-)zero input rather than NaN.
    Vec3 normalized() const {
        const Real len = length();
        return len > static_cast<Real>(1e-12) ? Vec3{x / len, y / len, z / len} : Vec3{};
    }
};

constexpr Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
constexpr Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
constexpr Vec3 operator*(Vec3 v, Real s) { return {v.x * s, v.y * s, v.z * s}; }
constexpr Vec3 operator*(Real s, Vec3 v) { return {v.x * s, v.y * s, v.z * s}; }
constexpr Vec3 operator/(Vec3 v, Real s) { return {v.x / s, v.y / s, v.z / s}; }

constexpr Real dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
// Right-handed cross product: x cross y = z. |a x b| is the area of the parallelogram they span, and
// a x b is perpendicular to both. Unlike 2D, this is a vector, because in 3D a rotation has an axis.
constexpr Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline Real distance(Vec3 a, Vec3 b) { return (a - b).length(); }

// 3x3 matrix, stored as three columns.
struct Mat3 {
    Vec3 cx{1, 0, 0}, cy{0, 1, 0}, cz{0, 0, 1};  // identity by default

    constexpr Mat3() = default;
    constexpr Mat3(Vec3 c0, Vec3 c1, Vec3 c2) : cx(c0), cy(c1), cz(c2) {}

    static constexpr Mat3 identity() { return {}; }
    static constexpr Mat3 diagonal(Real a, Real b, Real c) { return {{a, 0, 0}, {0, b, 0}, {0, 0, c}}; }
    // The matrix that does a cross product: skew(a) * b == cross(a, b).
    static constexpr Mat3 skew(Vec3 v) { return {{0, v.z, -v.y}, {-v.z, 0, v.x}, {v.y, -v.x, 0}}; }
    // a b^T: column j is a * b[j].
    static constexpr Mat3 outer(Vec3 a, Vec3 b) { return {a * b.x, a * b.y, a * b.z}; }

    // Element in row r, column c (both 0..2).
    constexpr Real at(int r, int c) const { return (c == 0 ? cx : c == 1 ? cy : cz)[r]; }

    constexpr Mat3 transposed() const {
        return {{cx.x, cy.x, cz.x}, {cx.y, cy.y, cz.y}, {cx.z, cy.z, cz.z}};
    }
    // Scalar triple product of the columns: the (signed) volume they span.
    constexpr Real det() const { return dot(cx, cross(cy, cz)); }

    // Each ROW of the inverse is a cross product of two columns of the original, divided by the
    // determinant. Returns the zero matrix when singular rather than dividing by zero.
    constexpr Mat3 inverse() const {
        const Real d = det();
        if (d == 0) return {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        const Vec3 r0 = cross(cy, cz) / d, r1 = cross(cz, cx) / d, r2 = cross(cx, cy) / d;
        return {{r0.x, r1.x, r2.x}, {r0.y, r1.y, r2.y}, {r0.z, r1.z, r2.z}};
    }
};

constexpr Vec3 operator*(const Mat3& m, Vec3 v) { return m.cx * v.x + m.cy * v.y + m.cz * v.z; }
constexpr Mat3 operator*(const Mat3& a, const Mat3& b) { return {a * b.cx, a * b.cy, a * b.cz}; }
constexpr Mat3 operator+(const Mat3& a, const Mat3& b) { return {a.cx + b.cx, a.cy + b.cy, a.cz + b.cz}; }
constexpr Mat3 operator-(const Mat3& a, const Mat3& b) { return {a.cx - b.cx, a.cy - b.cy, a.cz - b.cz}; }
constexpr Mat3 operator*(const Mat3& m, Real s) { return {m.cx * s, m.cy * s, m.cz * s}; }
constexpr Mat3 operator*(Real s, const Mat3& m) { return m * s; }

// Rotation as a unit quaternion. Composing rotations is a quaternion product, rotating a vector is
// cheap, and unlike three Euler angles it has no gimbal lock; unlike a 3x3 matrix it drifts back to a
// valid rotation with a single normalize().
struct Quat {
    Real x = 0, y = 0, z = 0, w = 1;  // identity by default

    constexpr Quat() = default;
    constexpr Quat(Real x_, Real y_, Real z_, Real w_) : x(x_), y(y_), z(z_), w(w_) {}

    static Quat from_axis_angle(Vec3 axis, Real angle) {
        const Vec3 a = axis.normalized() * std::sin(angle * static_cast<Real>(0.5));
        return {a.x, a.y, a.z, std::cos(angle * static_cast<Real>(0.5))};
    }

    constexpr Vec3 vector_part() const { return {x, y, z}; }
    constexpr Real length_sq() const { return x * x + y * y + z * z + w * w; }
    Real length() const { return std::sqrt(length_sq()); }
    // Back to unit length; the identity for a zero quaternion.
    Quat normalized() const {
        const Real len = length();
        return len > static_cast<Real>(1e-12) ? Quat{x / len, y / len, z / len, w / len} : Quat{};
    }
};

// Hamilton product. a * b means "rotate by b first, then by a" (the same order as matrices), so
// rotate(a * b, v) == rotate(a, rotate(b, v)). It does not commute.
constexpr Quat operator*(const Quat& a, const Quat& b) {
    return {a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
            a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
constexpr Quat operator*(const Quat& q, Real s) { return {q.x * s, q.y * s, q.z * s, q.w * s}; }
constexpr Quat operator+(const Quat& a, const Quat& b) { return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w}; }

// For a UNIT quaternion the inverse rotation is the conjugate.
constexpr Quat conjugate(const Quat& q) { return {-q.x, -q.y, -q.z, q.w}; }

// Rotate v: with u the vector part, v' = v + 2 w (u x v) + 2 u x (u x v). Cheaper than building a matrix.
constexpr Vec3 rotate(const Quat& q, Vec3 v) {
    const Vec3 u = q.vector_part();
    const Vec3 t = cross(u, v) * static_cast<Real>(2);
    return v + t * q.w + cross(u, t);
}
constexpr Vec3 inv_rotate(const Quat& q, Vec3 v) { return rotate(conjugate(q), v); }

// The 3x3 rotation matrix of a unit quaternion (columns are where the x, y, z axes go).
constexpr Mat3 to_mat3(const Quat& q) {
    const Real xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
    const Real xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
    const Real wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
    return {{1 - 2 * (yy + zz), 2 * (xy + wz), 2 * (xz - wy)},
            {2 * (xy - wz), 1 - 2 * (xx + zz), 2 * (yz + wx)},
            {2 * (xz + wy), 2 * (yz - wx), 1 - 2 * (xx + yy)}};
}

// The inverse conversion, for a proper rotation matrix. Branches on whichever of w, x, y, z is largest
// so it never divides by something close to zero.
inline Quat from_mat3(const Mat3& m) {
    const Real m00 = m.at(0, 0), m11 = m.at(1, 1), m22 = m.at(2, 2);
    const Real trace = m00 + m11 + m22;
    Quat q;
    if (trace > 0) {
        const Real s = std::sqrt(trace + 1) * 2;  // s = 4 w
        q = {(m.at(2, 1) - m.at(1, 2)) / s, (m.at(0, 2) - m.at(2, 0)) / s, (m.at(1, 0) - m.at(0, 1)) / s, s / 4};
    } else if (m00 > m11 && m00 > m22) {
        const Real s = std::sqrt(1 + m00 - m11 - m22) * 2;  // s = 4 x
        q = {s / 4, (m.at(0, 1) + m.at(1, 0)) / s, (m.at(0, 2) + m.at(2, 0)) / s, (m.at(2, 1) - m.at(1, 2)) / s};
    } else if (m11 > m22) {
        const Real s = std::sqrt(1 + m11 - m00 - m22) * 2;  // s = 4 y
        q = {(m.at(0, 1) + m.at(1, 0)) / s, s / 4, (m.at(1, 2) + m.at(2, 1)) / s, (m.at(0, 2) - m.at(2, 0)) / s};
    } else {
        const Real s = std::sqrt(1 + m22 - m00 - m11) * 2;  // s = 4 z
        q = {(m.at(0, 2) + m.at(2, 0)) / s, (m.at(1, 2) + m.at(2, 1)) / s, s / 4, (m.at(1, 0) - m.at(0, 1)) / s};
    }
    return q.normalized();
}

// The angle (0..pi) of the rotation a unit quaternion represents. Uses |w| so q and -q agree.
inline Real rotation_angle(const Quat& q) {
    return 2 * std::atan2(q.vector_part().length(), std::fabs(q.w));
}

// Advance an orientation by angular velocity `omega` (world frame, rad/s) over dt. The rate of change of
// a rotation quaternion is  dq/dt = 1/2 * (omega, 0) * q ; take one Euler step, then renormalise, which
// turns the small error of that step into a slightly wrong rotation angle instead of a non-rotation.
inline Quat integrate_orientation(const Quat& q, Vec3 omega, Real dt) {
    const Quat spin{omega.x, omega.y, omega.z, 0};
    return (q + (spin * q) * (static_cast<Real>(0.5) * dt)).normalized();
}

// Rotate then translate, like the 2D Transform.
struct Transform {
    Vec3 p;
    Quat q;
};
inline Vec3 apply(const Transform& t, Vec3 local) { return rotate(t.q, local) + t.p; }
inline Vec3 apply_inv(const Transform& t, Vec3 world) { return inv_rotate(t.q, world - t.p); }

}  // namespace phys3d
