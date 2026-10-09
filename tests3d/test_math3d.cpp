#include "test.hpp"
#include <phys3d/math.hpp>

#include <cstdint>

using namespace phys3d;

namespace {

constexpr double kEps = 1e-4;

struct Lcg {
    std::uint32_t s = 1;
    Real next() {  // uniform in [0, 1)
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Vec3 vec(Real extent = 3) { return {range(-extent, extent), range(-extent, extent), range(-extent, extent)}; }
    Quat rotation() {  // a random unit quaternion
        Quat q;
        do q = Quat{range(-1, 1), range(-1, 1), range(-1, 1), range(-1, 1)};
        while (q.length() < 0.2f);
        return q.normalized();
    }
    Mat3 matrix() { return {vec(1), vec(1), vec(1)}; }
};

#define CHECK_VEC(a, b, eps)                         \
    do {                                             \
        const Vec3 va_ = (a), vb_ = (b);             \
        CHECK_NEAR(va_.x, vb_.x, eps);               \
        CHECK_NEAR(va_.y, vb_.y, eps);               \
        CHECK_NEAR(va_.z, vb_.z, eps);               \
    } while (0)

#define CHECK_MAT(a, b, eps)                         \
    do {                                             \
        const Mat3 ma_ = (a), mb_ = (b);             \
        CHECK_VEC(ma_.cx, mb_.cx, eps);              \
        CHECK_VEC(ma_.cy, mb_.cy, eps);              \
        CHECK_VEC(ma_.cz, mb_.cz, eps);              \
    } while (0)

}  // namespace

// ---- Vec3 ---------------------------------------------------------------------------------------

TEST(vec3_arithmetic_and_length) {
    const Vec3 a{1, 2, 3}, b{4, -5, 6};
    CHECK_VEC(a + b, Vec3(5, -3, 9), kEps);
    CHECK_VEC(a - b, Vec3(-3, 7, -3), kEps);
    CHECK_VEC(2.0f * a, Vec3(2, 4, 6), kEps);
    CHECK_VEC(b / 2.0f, Vec3(2, -2.5f, 3), kEps);
    CHECK_VEC(-a, Vec3(-1, -2, -3), kEps);
    CHECK_NEAR(dot(a, b), 4 - 10 + 18, kEps);
    CHECK_NEAR(Vec3(2, 3, 6).length(), 7, kEps);  // 4 + 9 + 36 = 49
    CHECK_NEAR(Vec3(2, 3, 6).normalized().length(), 1, kEps);
    CHECK(Vec3{}.normalized().length_sq() == 0);  // no NaN from the zero vector
    Vec3 v{1, 2, 3};
    v[1] = 9;
    CHECK(v.y == 9 && v[0] == 1 && v[2] == 3);
}

TEST(cross_product_is_right_handed) {
    CHECK_VEC(cross({1, 0, 0}, {0, 1, 0}), Vec3(0, 0, 1), kEps);  // x cross y = z
    CHECK_VEC(cross({0, 1, 0}, {0, 0, 1}), Vec3(1, 0, 0), kEps);  // y cross z = x
    CHECK_VEC(cross({0, 0, 1}, {1, 0, 0}), Vec3(0, 1, 0), kEps);  // z cross x = y
    CHECK_VEC(cross({0, 1, 0}, {1, 0, 0}), Vec3(0, 0, -1), kEps);
    CHECK_VEC(cross({1, 2, 3}, {4, 5, 6}), Vec3(-3, 6, -3), kEps);  // a textbook example
}

TEST(cross_product_identities_hold_for_random_vectors) {
    Lcg rng;
    for (int i = 0; i < 300; ++i) {
        const Vec3 a = rng.vec(), b = rng.vec(), c = rng.vec();
        CHECK_VEC(cross(a, b), -cross(b, a), 1e-4);                      // anticommutative
        CHECK_NEAR(dot(cross(a, b), a), 0, 2e-4);                        // perpendicular to both inputs
        CHECK_NEAR(dot(cross(a, b), b), 0, 2e-4);
        // Lagrange: |a x b|^2 = |a|^2 |b|^2 - (a.b)^2
        CHECK_NEAR(cross(a, b).length_sq(), a.length_sq() * b.length_sq() - dot(a, b) * dot(a, b), 2e-2);
        // the "BAC-CAB" rule: a x (b x c) = b (a.c) - c (a.b)
        CHECK_VEC(cross(a, cross(b, c)), b * dot(a, c) - c * dot(a, b), 2e-3);
        // the scalar triple product is unchanged by cycling its arguments
        CHECK_NEAR(dot(a, cross(b, c)), dot(b, cross(c, a)), 2e-3);
        // Jacobi identity
        CHECK_VEC(cross(a, cross(b, c)) + cross(b, cross(c, a)) + cross(c, cross(a, b)), Vec3(), 5e-3);
    }
}

// ---- Mat3 ---------------------------------------------------------------------------------------

TEST(mat3_known_inverse) {
    // Rows [1 0 5; 2 1 6; 3 4 0]: the classic determinant-1 example, inverse [-24 20 -5; 18 -15 4; 5 -4 1].
    const Mat3 m{{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};  // columns
    CHECK_NEAR(m.det(), 1, kEps);
    const Mat3 inv = m.inverse();
    const Real expected[3][3] = {{-24, 20, -5}, {18, -15, 4}, {5, -4, 1}};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) CHECK_NEAR(inv.at(r, c), expected[r][c], 1e-3);
    CHECK_MAT(m * inv, Mat3::identity(), 1e-3);
}

TEST(mat3_layout_matches_the_documentation) {
    const Mat3 m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};  // columns
    CHECK(m.at(0, 0) == 1 && m.at(1, 0) == 2 && m.at(2, 0) == 3);  // first column
    CHECK(m.at(0, 1) == 4 && m.at(0, 2) == 7);                     // first row is (1, 4, 7)
    CHECK_VEC(m * Vec3(1, 0, 0), Vec3(1, 2, 3), kEps);             // M e_x is the first column
    CHECK_VEC(m * Vec3(1, 1, 1), Vec3(12, 15, 18), kEps);
    CHECK_VEC(m.transposed().cx, Vec3(1, 4, 7), kEps);
}

TEST(mat3_algebra_for_random_matrices) {
    Lcg rng;
    int checked = 0;
    for (int i = 0; i < 300; ++i) {
        const Mat3 a = rng.matrix(), b = rng.matrix();
        CHECK_NEAR((a * b).det(), a.det() * b.det(), 5e-3);                       // det of a product
        CHECK_MAT((a * b).transposed(), b.transposed() * a.transposed(), 1e-4);   // transpose reverses order
        CHECK_MAT((a + b) * 2.0f, a * 2.0f + b * 2.0f, 1e-4);
        const Vec3 v = rng.vec();
        CHECK_VEC((a * b) * v, a * (b * v), 1e-3);                                // associativity on a vector
        if (std::fabs(a.det()) > 0.15f) {                                         // well-conditioned only
            CHECK_MAT(a * a.inverse(), Mat3::identity(), 2e-3);
            CHECK_MAT(a.inverse() * a, Mat3::identity(), 2e-3);
            ++checked;
        }
    }
    CHECK(checked > 100);
}

TEST(singular_matrix_inverts_to_zero_not_nan) {
    const Mat3 flat{{1, 2, 3}, {2, 4, 6}, {0, 1, 0}};  // second column is twice the first
    const Mat3 z = flat.inverse();
    CHECK(z.cx.length_sq() == 0 && z.cy.length_sq() == 0 && z.cz.length_sq() == 0);
}

TEST(skew_matrix_does_a_cross_product) {
    Lcg rng;
    for (int i = 0; i < 100; ++i) {
        const Vec3 a = rng.vec(), b = rng.vec();
        const Mat3 s = Mat3::skew(a);
        CHECK_VEC(s * b, cross(a, b), 1e-4);
        CHECK_MAT(s.transposed(), s * -1.0f, 1e-5);  // antisymmetric
        // [a]x^2 = a a^T - |a|^2 I
        CHECK_MAT(s * s, Mat3::outer(a, a) - Mat3::identity() * a.length_sq(), 1e-3);
    }
}

// ---- Quat ---------------------------------------------------------------------------------------

TEST(quaternion_rotation_follows_the_right_hand_rule) {
    const Real q90 = kPi / 2;
    CHECK_VEC(rotate(Quat::from_axis_angle({0, 0, 1}, q90), Vec3(1, 0, 0)), Vec3(0, 1, 0), kEps);  // about z: x -> y
    CHECK_VEC(rotate(Quat::from_axis_angle({1, 0, 0}, q90), Vec3(0, 1, 0)), Vec3(0, 0, 1), kEps);  // about x: y -> z
    CHECK_VEC(rotate(Quat::from_axis_angle({0, 1, 0}, q90), Vec3(0, 0, 1)), Vec3(1, 0, 0), kEps);  // about y: z -> x
    CHECK_VEC(rotate(Quat::from_axis_angle({0, 0, 1}, kPi), Vec3(1, 2, 3)), Vec3(-1, -2, 3), kEps);  // half turn
    // The rotation axis is fixed, and a non-unit axis is normalised for us.
    CHECK_VEC(rotate(Quat::from_axis_angle({3, 0, 4}, 1.1f), Vec3(3, 0, 4)), Vec3(3, 0, 4), 1e-4);
    CHECK_VEC(rotate(Quat{}, Vec3(1, 2, 3)), Vec3(1, 2, 3), kEps);  // identity
}

TEST(quaternion_composition_order) {
    const Quat about_x = Quat::from_axis_angle({1, 0, 0}, kPi / 2), about_y = Quat::from_axis_angle({0, 1, 0}, kPi / 2);
    const Vec3 v{0, 0, 1};
    // a * b = "b first, then a".
    CHECK_VEC(rotate(about_y * about_x, v), rotate(about_y, rotate(about_x, v)), kEps);
    // Rotations do not commute: x then y is not y then x.
    const Vec3 xy = rotate(about_y * about_x, Vec3(0, 1, 0)), yx = rotate(about_x * about_y, Vec3(0, 1, 0));
    CHECK((xy - yx).length() > 1.0f);
}

TEST(quaternion_properties_for_random_rotations) {
    Lcg rng;
    for (int i = 0; i < 300; ++i) {
        const Quat a = rng.rotation(), b = rng.rotation(), c = rng.rotation();
        const Vec3 v = rng.vec();
        CHECK_NEAR(rotate(a, v).length(), v.length(), 1e-4);                         // rotations preserve length
        CHECK_VEC(rotate(a * b, v), rotate(a, rotate(b, v)), 2e-4);                  // composition
        CHECK_VEC(inv_rotate(a, rotate(a, v)), v, 2e-4);                             // inverse
        CHECK_VEC(rotate(Quat{-a.x, -a.y, -a.z, -a.w}, v), rotate(a, v), 1e-4);      // q and -q are one rotation
        CHECK_VEC(rotate((a * b) * c, v), rotate(a * (b * c), v), 2e-4);             // associative
        CHECK_NEAR((a * b).length(), 1, 1e-4);                                       // products of unit quaternions stay unit
        const Mat3 m = to_mat3(a);
        CHECK_VEC(m * v, rotate(a, v), 2e-4);                                        // the matrix does the same thing
        CHECK_MAT(m.transposed() * m, Mat3::identity(), 1e-4);                       // orthonormal...
        CHECK_NEAR(m.det(), 1, 1e-4);                                                // ...and not a reflection
        CHECK_MAT(to_mat3(a * b), to_mat3(a) * to_mat3(b), 2e-4);                    // matrices compose the same way
    }
}

TEST(quaternion_matrix_round_trip_covers_every_branch) {
    Lcg rng;
    int branches[4] = {0, 0, 0, 0};
    for (int i = 0; i < 400; ++i) {
        const Quat q = rng.rotation();
        const Mat3 m = to_mat3(q);
        const Real tr = m.at(0, 0) + m.at(1, 1) + m.at(2, 2);
        ++branches[tr > 0 ? 0 : (m.at(0, 0) > m.at(1, 1) && m.at(0, 0) > m.at(2, 2)) ? 1 : (m.at(1, 1) > m.at(2, 2)) ? 2 : 3];
        const Quat back = from_mat3(m);
        const Vec3 v = rng.vec();
        CHECK_VEC(rotate(back, v), rotate(q, v), 3e-4);  // the same rotation (maybe the negated quaternion)
    }
    for (int b : branches) CHECK(b > 5);  // the test really exercised all four cases

    // Half turns about each axis land in the awkward branches.
    for (int axis = 0; axis < 3; ++axis) {
        Vec3 ax;
        ax[axis] = 1;
        const Quat q = Quat::from_axis_angle(ax, kPi);
        CHECK_VEC(rotate(from_mat3(to_mat3(q)), Vec3(1, 2, 3)), rotate(q, Vec3(1, 2, 3)), 1e-4);
    }
}

TEST(rotation_angle_is_recovered) {
    Lcg rng;
    for (int i = 0; i < 100; ++i) {
        const Real angle = rng.range(0.05f, 3.1f);
        const Quat q = Quat::from_axis_angle(rng.vec(), angle);
        CHECK_NEAR(rotation_angle(q), angle, 1e-3);
        CHECK_NEAR(rotation_angle(Quat{-q.x, -q.y, -q.z, -q.w}), angle, 1e-3);  // same for -q
    }
    CHECK_NEAR(rotation_angle(Quat{}), 0, kEps);
}

// Constant angular velocity: after one second the body has turned by exactly |omega| radians about omega.
TEST(orientation_integration_turns_at_the_angular_velocity) {
    const Vec3 omega{1.2f, -2.0f, 2.4f};  // |omega| = 3.4
    Quat q;
    const Real dt = 1.0f / 120.0f;
    for (int i = 0; i < 120; ++i) q = integrate_orientation(q, omega, dt);
    const Quat exact = Quat::from_axis_angle(omega, omega.length());
    const Vec3 probe{1, 0.5f, -0.3f};
    CHECK_VEC(rotate(q, probe), rotate(exact, probe), 2e-3);
    CHECK_NEAR(q.length(), 1, 1e-5);
}

TEST(orientation_integration_stays_a_valid_rotation_over_a_long_run) {
    Lcg rng;
    Quat q;
    for (int i = 0; i < 20000; ++i) q = integrate_orientation(q, Vec3{3, 1, -2} + rng.vec(0.5f), 1.0f / 60.0f);
    CHECK_NEAR(q.length(), 1, 1e-5);
    CHECK_NEAR(to_mat3(q).det(), 1, 1e-4);
}

TEST(transform_round_trip_and_rigidity) {
    Lcg rng;
    for (int i = 0; i < 100; ++i) {
        const Transform t{rng.vec(5), rng.rotation()};
        const Vec3 a = rng.vec(), b = rng.vec();
        CHECK_VEC(apply_inv(t, apply(t, a)), a, 3e-4);
        CHECK_NEAR(distance(apply(t, a), apply(t, b)), distance(a, b), 3e-4);  // rigid: distances preserved
    }
}
