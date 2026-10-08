#include "test.hpp"
#include <phys/math.hpp>

using namespace phys;
constexpr double kEps = 1e-5;

TEST(vec_arithmetic) {
    Vec2 a{1, 2}, b{3, -4};
    Vec2 s = a + b, d = a - b, m = 2.0f * a;
    CHECK_NEAR(s.x, 4, kEps); CHECK_NEAR(s.y, -2, kEps);
    CHECK_NEAR(d.x, -2, kEps); CHECK_NEAR(d.y, 6, kEps);
    CHECK_NEAR(m.x, 2, kEps); CHECK_NEAR(m.y, 4, kEps);
    CHECK_NEAR(dot(a, b), -5, kEps);
}

TEST(vec_length_and_normalize) {
    Vec2 v{3, 4};
    CHECK_NEAR(v.length(), 5, kEps);
    CHECK_NEAR(v.normalized().length(), 1, kEps);
    Vec2 z = Vec2{}.normalized();  // must not produce NaN
    CHECK(z.x == 0 && z.y == 0);
}

TEST(cross_products) {
    // x cross y = +1 (counter-clockwise)
    CHECK_NEAR(cross(Vec2{1, 0}, Vec2{0, 1}), 1, kEps);
    CHECK_NEAR(cross(Vec2{0, 1}, Vec2{1, 0}), -1, kEps);
    // w x r: spinning CCW at offset (1,0) moves the point toward +y
    Vec2 v = cross(2.0f, Vec2{1, 0});
    CHECK_NEAR(v.x, 0, kEps); CHECK_NEAR(v.y, 2, kEps);
    // r x w is the negative of w x r
    Vec2 u = cross(Vec2{1, 0}, 2.0f);
    CHECK_NEAR(u.x, -v.x, kEps); CHECK_NEAR(u.y, -v.y, kEps);
}

TEST(rotation) {
    Rot q(kPi / 2);
    Vec2 r = rotate(q, Vec2{1, 0});
    CHECK_NEAR(r.x, 0, kEps); CHECK_NEAR(r.y, 1, kEps);
    CHECK_NEAR(q.angle(), kPi / 2, kEps);
    Vec2 back = inv_rotate(q, r);
    CHECK_NEAR(back.x, 1, kEps); CHECK_NEAR(back.y, 0, kEps);
}

TEST(mat2_inverse) {
    Mat2 m{{2, 1}, {1, 3}};  // columns
    Mat2 inv = m.inverse();
    Vec2 v{5, -7};
    Vec2 round_trip = inv * (m * v);
    CHECK_NEAR(round_trip.x, 5, kEps); CHECK_NEAR(round_trip.y, -7, kEps);
    CHECK_NEAR(m.det(), 5, kEps);
    Mat2 singular{{1, 2}, {2, 4}};
    Mat2 z = singular.inverse();  // zero matrix, not NaN
    CHECK(z.ex.x == 0 && z.ey.y == 0);
}

TEST(transform_round_trip) {
    Transform t(Vec2{3, -2}, 0.7f);
    Vec2 local{1.5f, 0.25f};
    Vec2 back = apply_inv(t, apply(t, local));
    CHECK_NEAR(back.x, local.x, kEps); CHECK_NEAR(back.y, local.y, kEps);
    // Rigid: distances are preserved
    Vec2 a = apply(t, Vec2{0, 0}), b = apply(t, Vec2{3, 4});
    CHECK_NEAR(distance(a, b), 5, kEps);
}
