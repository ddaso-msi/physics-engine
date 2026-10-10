#include "test.hpp"
#include <phys3d/inertia.hpp>

#include <cstdint>
#include <utility>

using namespace phys3d;

namespace {

struct Lcg {
    std::uint32_t s = 5;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Quat rotation() {
        Quat q;
        do q = Quat{range(-1, 1), range(-1, 1), range(-1, 1), range(-1, 1)};
        while (q.length() < 0.2f);
        return q.normalized();
    }
};

#define CHECK_MAT(a, b, eps)                                  \
    do {                                                      \
        const Mat3 ma_ = (a), mb_ = (b);                      \
        for (int r_ = 0; r_ < 3; ++r_)                        \
            for (int c_ = 0; c_ < 3; ++c_)                    \
                CHECK_NEAR(ma_.at(r_, c_), mb_.at(r_, c_), eps); \
    } while (0)

// Inertia tensor about `pivot` of a uniform body, by brute-force integration: the body fills the points
// of an n^3 midpoint grid over `lo..hi` for which inside(p) is true, each point carrying an equal share of
// the mass. I_ij = sum m_k (|r|^2 delta_ij - r_i r_j), with r measured from the pivot. An independent
// check on the formulas that does not use them.
template <class Inside>
Mat3 integrate_inertia(Real mass, Vec3 lo, Vec3 hi, Vec3 pivot, int n, Inside inside) {
    int count = 0;
    Mat3 sum = Mat3::diagonal(0, 0, 0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                const Vec3 p{lo.x + (hi.x - lo.x) * (static_cast<Real>(i) + 0.5f) / static_cast<Real>(n),
                             lo.y + (hi.y - lo.y) * (static_cast<Real>(j) + 0.5f) / static_cast<Real>(n),
                             lo.z + (hi.z - lo.z) * (static_cast<Real>(k) + 0.5f) / static_cast<Real>(n)};
                if (!inside(p)) continue;
                const Vec3 r = p - pivot;
                sum = sum + (Mat3::identity() * r.length_sq() - Mat3::outer(r, r));
                ++count;
            }
    return sum * (mass / static_cast<Real>(count));
}

}  // namespace

TEST(sphere_and_box_inertia_match_the_textbook_formulas) {
    const Mat3 s = sphere_inertia(5, 2);
    CHECK_NEAR(s.at(0, 0), 0.4 * 5 * 4, 1e-4);  // 2/5 m r^2
    CHECK_NEAR(s.at(1, 1), 8, 1e-4);
    CHECK_NEAR(s.at(2, 2), 8, 1e-4);
    CHECK_NEAR(s.at(0, 1), 0, 1e-6);

    // A 2 x 4 x 6 box (half extents 1, 2, 3), mass 12: Ixx = m (4^2 + 6^2)/12 = 52, and so on.
    const Mat3 b = box_inertia(12, {1, 2, 3});
    CHECK_NEAR(b.at(0, 0), 12.0 * (16 + 36) / 12, 1e-3);
    CHECK_NEAR(b.at(1, 1), 12.0 * (4 + 36) / 12, 1e-3);
    CHECK_NEAR(b.at(2, 2), 12.0 * (4 + 16) / 12, 1e-3);
    CHECK_NEAR(b.at(0, 1) + b.at(0, 2) + b.at(1, 2), 0, 1e-6);  // principal axes: no off-diagonal terms
}

TEST(box_inertia_obeys_the_triangle_inequality_of_a_real_body) {
    // For any physical body each principal moment is at most the sum of the other two.
    Lcg rng;
    for (int i = 0; i < 100; ++i) {
        const Mat3 b = box_inertia(rng.range(0.5f, 20), {rng.range(0.1f, 3), rng.range(0.1f, 3), rng.range(0.1f, 3)});
        const Real x = b.at(0, 0), y = b.at(1, 1), z = b.at(2, 2);
        CHECK(x <= y + z + 1e-4f && y <= x + z + 1e-4f && z <= x + y + 1e-4f);
    }
}

TEST(inertia_formulas_agree_with_numerical_integration) {
    const Vec3 half{0.8f, 0.5f, 1.3f};
    const Mat3 box_numeric = integrate_inertia(3, -half, half, {0, 0, 0}, 40, [](Vec3) { return true; });
    CHECK_MAT(box_numeric, box_inertia(3, half), 5e-3);

    const Real r = 1.1f;
    const Mat3 sphere_numeric = integrate_inertia(3, {-r, -r, -r}, {r, r, r}, {0, 0, 0}, 60,
                                                  [r](Vec3 p) { return p.length_sq() <= r * r; });
    const Mat3 sphere_exact = sphere_inertia(3, r);
    CHECK_NEAR(sphere_numeric.at(0, 0), sphere_exact.at(0, 0), 0.02 * sphere_exact.at(0, 0));  // a few % from the staircase surface
    CHECK_NEAR(sphere_numeric.at(1, 1), sphere_exact.at(1, 1), 0.02 * sphere_exact.at(1, 1));
    CHECK_NEAR(sphere_numeric.at(0, 1), 0, 0.01);
}

TEST(capsule_inertia_agrees_with_numerical_integration) {
    for (const auto& [half_length, radius] : {std::pair<Real, Real>{0.8f, 0.3f}, {0.2f, 0.5f}, {1.5f, 0.1f}}) {
        const Real mass = 3, reach = half_length + radius;
        auto inside = [&](Vec3 p) {
            const Real y = p.y < -half_length ? -half_length : p.y > half_length ? half_length : p.y;
            return (p - Vec3{0, y, 0}).length_sq() <= radius * radius;
        };
        const Mat3 slow = integrate_inertia(mass, {-radius, -reach, -radius}, {radius, reach, radius}, {0, 0, 0}, 120, inside);
        CHECK_MAT(capsule_inertia(mass, half_length, radius), slow, 0.004 * slow.at(0, 0));
        // And the volume, the same way: the share of the grid's box that is inside.
        int in = 0;
        const int n = 100;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k)
                    in += inside({radius * (2 * (i + 0.5f) / n - 1), reach * (2 * (j + 0.5f) / n - 1), radius * (2 * (k + 0.5f) / n - 1)});
        const double box_volume = 8.0 * radius * reach * radius;
        CHECK_NEAR(capsule_volume(half_length, radius), box_volume * in / (double(n) * n * n), 0.005 * box_volume);
    }
}

TEST(capsule_inertia_has_the_right_limits) {
    // No cylinder at all: a sphere.
    CHECK_MAT(capsule_inertia(2, 0, 0.7f), sphere_inertia(2, 0.7f), 1e-5);
    // Very thin: a rod of length 2h, m L^2 / 12 = m h^2 / 3 across and nothing about its axis.
    const Mat3 rod = capsule_inertia(2, 1.5f, 1e-3f);
    CHECK_NEAR(rod.at(0, 0), 2 * 1.5 * 1.5 / 3, 3e-3);
    CHECK_NEAR(rod.at(2, 2), rod.at(0, 0), 1e-6);
    CHECK(rod.at(1, 1) < 1e-5f);
    // Harder to turn end over end than to spin about its axis.
    const Mat3 c = capsule_inertia(1, 0.5f, 0.25f);
    CHECK(c.at(0, 0) > c.at(1, 1));
}

TEST(parallel_axis_theorem) {
    // A box pivoted about one of its corners. The offset from the centre of mass to the pivot is d = (a, b, c)
    // (the half extents), and the tensor gains off-diagonal terms, so compare with the brute-force integral
    // about that corner instead of a closed form.
    const Vec3 half{0.8f, 0.5f, 1.3f};
    const Vec3 corner = half;
    const Mat3 numeric = integrate_inertia(3, -half, half, corner, 40, [](Vec3) { return true; });
    const Mat3 shifted = parallel_axis(box_inertia(3, half), 3, corner);
    CHECK_MAT(shifted, numeric, 2e-2);
    // Off-diagonal terms appear: I_xy = -m d_x d_y
    CHECK_NEAR(shifted.at(0, 1), -3.0 * 0.8 * 0.5, 1e-3);
    CHECK_NEAR(shifted.at(0, 1), shifted.at(1, 0), 1e-6);  // symmetric
    // A zero offset changes nothing.
    CHECK_MAT(parallel_axis(box_inertia(3, half), 3, {0, 0, 0}), box_inertia(3, half), 1e-6);
}

TEST(rotating_a_box_inertia_into_the_world_frame) {
    const Vec3 half{1, 2, 3};
    const Mat3 body = box_inertia(12, half);

    // A quarter turn about z swaps the x and y axes: Ixx and Iyy trade places.
    const Mat3 world = to_world_frame(to_mat3(Quat::from_axis_angle({0, 0, 1}, kPi / 2)), body);
    CHECK_NEAR(world.at(0, 0), body.at(1, 1), 1e-3);
    CHECK_NEAR(world.at(1, 1), body.at(0, 0), 1e-3);
    CHECK_NEAR(world.at(2, 2), body.at(2, 2), 1e-3);

    // For any rotation: still symmetric, and trace, determinant (the "size" of the tensor) are unchanged.
    Lcg rng;
    for (int i = 0; i < 100; ++i) {
        const Mat3 r = to_mat3(rng.rotation());
        const Mat3 w = to_world_frame(r, body);
        CHECK_MAT(w, w.transposed(), 1e-3);
        CHECK_NEAR(w.at(0, 0) + w.at(1, 1) + w.at(2, 2), body.at(0, 0) + body.at(1, 1) + body.at(2, 2), 1e-3);
        CHECK_NEAR(w.det(), body.det(), 1e-4 * body.det());  // relative: the determinant here is about 40000
        // The inverse of the rotated tensor is the rotated inverse: what the solver actually stores.
        CHECK_MAT(to_world_frame(r, body.inverse()), w.inverse(), 1e-3);
    }
}

TEST(a_cube_has_the_same_inertia_whichever_way_it_is_turned) {
    const Mat3 cube = box_inertia(5, {0.7f, 0.7f, 0.7f});
    Lcg rng;
    for (int i = 0; i < 100; ++i) CHECK_MAT(to_world_frame(to_mat3(rng.rotation()), cube), cube, 1e-3);
    // A sphere likewise.
    const Mat3 ball = sphere_inertia(5, 0.9f);
    for (int i = 0; i < 100; ++i) CHECK_MAT(to_world_frame(to_mat3(rng.rotation()), ball), ball, 1e-3);
}

// L = I w is not parallel to w unless w lies along a principal axis: the signature of 3D rotation, and the
// reason the solver cannot treat angular velocity as a single number like the 2D engine does.
TEST(angular_momentum_is_not_parallel_to_angular_velocity_in_general) {
    const Mat3 box = box_inertia(12, {1, 2, 3});
    const Vec3 along_axis{0, 0, 2};
    const Vec3 tilted{1, 1, 1};
    const Vec3 l_axis = box * along_axis, l_tilted = box * tilted;
    CHECK_NEAR(cross(l_axis, along_axis).length(), 0, 1e-4);   // parallel
    CHECK(cross(l_tilted, tilted).length() > 1.0f);            // not parallel
}
