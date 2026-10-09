#include "test.hpp"
#include <phys3d/collision.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

using namespace phys3d;

namespace {

constexpr double kEps = 1e-3;

struct Lcg {
    std::uint32_t s = 21;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Vec3 vec(Real extent) { return {range(-extent, extent), range(-extent, extent), range(-extent, extent)}; }
    Quat rotation() {
        Quat q;
        do q = Quat{range(-1, 1), range(-1, 1), range(-1, 1), range(-1, 1)};
        while (q.length() < 0.2f);
        return q.normalized();
    }
};

#define CHECK_VEC(a, b, eps)                         \
    do {                                             \
        const Vec3 va_ = (a), vb_ = (b);             \
        CHECK_NEAR(va_.x, vb_.x, eps);               \
        CHECK_NEAR(va_.y, vb_.y, eps);               \
        CHECK_NEAR(va_.z, vb_.z, eps);               \
    } while (0)

Body box(Vec3 half, Vec3 pos, Quat q = {}) { return Body::solid_box(half, 1, pos, q); }
Body ball(Real r, Vec3 pos) { return Body::solid_sphere(r, 1, pos); }

Body random_body(Lcg& rng) {
    const Vec3 pos = rng.vec(1.4f);
    if (rng.next() < 0.3f) return ball(rng.range(0.3f, 1.0f), pos);
    return box({rng.range(0.3f, 1.2f), rng.range(0.3f, 1.2f), rng.range(0.3f, 1.2f)}, pos, rng.rotation());
}

// Is p inside the body's shape grown by `margin` on every side?
bool contains_inflated(const Body& b, Vec3 p, Real margin) {
    const Vec3 local = inv_rotate(b.q, p - b.pos);
    if (b.shape.type == Shape::Type::Sphere) return local.length() <= b.shape.radius + margin;
    const Vec3 h = b.shape.half_extents;
    return std::fabs(local.x) <= h.x + margin && std::fabs(local.y) <= h.y + margin && std::fabs(local.z) <= h.z + margin;
}

// A uniformly random point inside the body's shape.
Vec3 point_inside(const Body& b, Lcg& rng) {
    if (b.shape.type == Shape::Type::Sphere) {
        Vec3 p;
        do p = rng.vec(1);
        while (p.length_sq() > 1);
        return b.pos + p * b.shape.radius;
    }
    const Vec3 h = b.shape.half_extents;
    return apply(b.transform(), {rng.range(-h.x, h.x), rng.range(-h.y, h.y), rng.range(-h.z, h.z)});
}

std::vector<std::uint32_t> sorted_ids(const Manifold& m) {
    std::vector<std::uint32_t> ids;
    for (int i = 0; i < m.count; ++i) ids.push_back(m.points[i].id);
    std::sort(ids.begin(), ids.end());
    return ids;
}

// Area of the polygon formed by a manifold's points, seen along y (they must lie in a horizontal plane).
double patch_area_xz(const Manifold& m) {
    Vec3 centre;
    for (int i = 0; i < m.count; ++i) centre += m.points[i].point;
    centre = centre / static_cast<Real>(m.count);
    std::vector<Vec3> p;
    for (int i = 0; i < m.count; ++i) p.push_back(m.points[i].point - centre);
    std::sort(p.begin(), p.end(), [](Vec3 a, Vec3 b) { return std::atan2(a.z, a.x) < std::atan2(b.z, b.x); });
    double area = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        const Vec3 a = p[i], b = p[(i + 1) % p.size()];
        area += 0.5 * double(a.x * b.z - a.z * b.x);
    }
    return std::fabs(area);
}

}  // namespace

// ---- spheres --------------------------------------------------------------------------------------

TEST(spheres_overlap_and_separate) {
    Manifold m;
    CHECK(collide(ball(1.0f, {0, 0, 0}), ball(0.5f, {0, 1.2f, 0}), m));
    CHECK_VEC(m.normal, Vec3(0, 1, 0), kEps);
    CHECK_NEAR(m.depth, 0.3, kEps);
    CHECK(m.count == 1);
    CHECK_VEC(m.points[0].point, Vec3(0, 0.85f, 0), kEps);  // halfway through the overlap [0.7, 1.0]
    CHECK(!collide(ball(1.0f, {0, 0, 0}), ball(0.5f, {0, 1.6f, 0}), m));

    // Along a diagonal, to make sure nothing is special about the axes.
    CHECK(collide(ball(1.0f, {1, 1, 1}), ball(1.0f, {2, 2, 2}), m));
    CHECK_VEC(m.normal, Vec3(1, 1, 1).normalized(), kEps);
    CHECK_NEAR(m.depth, 2 - std::sqrt(3.0), kEps);

    CHECK(collide(ball(1.0f, {2, 3, 4}), ball(0.5f, {2, 3, 4}), m));  // concentric: no NaN
    CHECK_NEAR(m.normal.length(), 1, kEps);
    CHECK_NEAR(m.depth, 1.5, kEps);
}

// ---- box vs sphere: one case per feature of the box ---------------------------------------------------

TEST(sphere_against_a_box_face_edge_and_corner) {
    const Body cube = box({1, 1, 1}, {0, 0, 0});
    Manifold m;

    CHECK(collide(cube, ball(1.0f, {1.8f, 0.3f, -0.2f}), m));  // face
    CHECK_VEC(m.normal, Vec3(1, 0, 0), kEps);
    CHECK_NEAR(m.depth, 0.2, kEps);
    CHECK_VEC(m.points[0].point, Vec3(0.9f, 0.3f, -0.2f), kEps);  // between the face (x=1) and the sphere's skin (x=0.8)

    CHECK(collide(cube, ball(1.0f, {1.5f, 1.5f, 0}), m));  // edge: nearest point of the box is (1,1,0)
    CHECK_VEC(m.normal, Vec3(1, 1, 0).normalized(), kEps);
    CHECK_NEAR(m.depth, 1 - std::sqrt(0.5), kEps);

    CHECK(collide(cube, ball(1.0f, {1.5f, 1.5f, 1.5f}), m));  // corner
    CHECK_VEC(m.normal, Vec3(1, 1, 1).normalized(), kEps);
    CHECK_NEAR(m.depth, 1 - std::sqrt(0.75), kEps);

    // Within a radius of three faces at once, yet clear of the corner itself.
    CHECK(!collide(cube, ball(1.0f, {1.7f, 1.7f, 1.7f}), m));

    const Manifold face_m = [&] { Manifold x; collide(cube, ball(1.0f, {1.8f, 0, 0}), x); return x; }();
    const Manifold edge_m = [&] { Manifold x; collide(cube, ball(1.0f, {1.5f, 1.5f, 0}), x); return x; }();
    CHECK(face_m.points[0].id != edge_m.points[0].id);  // different features, different ids
}

TEST(sphere_centre_inside_a_box_leaves_by_the_nearest_face) {
    Manifold m;
    CHECK(collide(box({1, 2, 3}, {0, 0, 0}), ball(0.5f, {0.7f, 0.1f, -0.2f}), m));
    CHECK_VEC(m.normal, Vec3(1, 0, 0), kEps);  // 0.3 from the +x face, further from every other
    CHECK_NEAR(m.depth, 0.3 + 0.5, kEps);
    CHECK(collide(box({1, 2, 3}, {0, 0, 0}), ball(0.5f, {0.1f, -1.9f, 0.2f}), m));
    CHECK_VEC(m.normal, Vec3(0, -1, 0), kEps);
}

TEST(box_sphere_respects_the_box_orientation_and_argument_order) {
    // A 4 x 2 x 1 box turned a quarter turn about z: its long axis now points up.
    const Body tall = box({2, 1, 0.5f}, {0, 0, 0}, Quat::from_axis_angle({0, 0, 1}, kPi / 2));
    const Body marble = ball(0.5f, {0, 2.3f, 0});
    Manifold bs, sb;
    CHECK(collide(tall, marble, bs));
    CHECK_VEC(bs.normal, Vec3(0, 1, 0), kEps);
    CHECK_NEAR(bs.depth, 0.2, kEps);
    CHECK(collide(marble, tall, sb));  // sphere first: same contact, normal reversed
    CHECK_VEC(sb.normal, Vec3(0, -1, 0), kEps);
    CHECK_NEAR(sb.depth, bs.depth, kEps);
    CHECK_VEC(sb.points[0].point, bs.points[0].point, kEps);
}

// ---- box vs box ---------------------------------------------------------------------------------------

TEST(box_resting_on_a_box_gives_a_four_point_patch) {
    // A small crate on a big one, sunk 0.1 into its top. The contact patch is the small crate's whole
    // bottom face: four points at its corners, each halfway through the overlap.
    Manifold m;
    CHECK(collide(box({1, 1, 1}, {0, 0, 0}), box({0.5f, 0.5f, 0.5f}, {0.2f, 1.4f, -0.1f}), m));
    CHECK_VEC(m.normal, Vec3(0, 1, 0), kEps);
    CHECK_NEAR(m.depth, 0.1, kEps);
    CHECK(m.count == 4);
    for (int i = 0; i < m.count; ++i) {
        CHECK_NEAR(m.points[i].point.y, 0.95, kEps);
        CHECK_NEAR(m.points[i].depth, 0.1, kEps);
        CHECK_NEAR(std::fabs(m.points[i].point.x - 0.2f), 0.5, kEps);
        CHECK_NEAR(std::fabs(m.points[i].point.z + 0.1f), 0.5, kEps);
    }
    CHECK_NEAR(patch_area_xz(m), 1.0, 1e-2);
    const auto ids = sorted_ids(m);
    CHECK(std::adjacent_find(ids.begin(), ids.end()) == ids.end());  // all different
}

TEST(overhanging_box_is_clipped_to_the_face_below) {
    // The small crate hangs 0.3 over the edge: its bottom face spans x 0.3..1.3, the face below ends at 1.
    Manifold m;
    CHECK(collide(box({1, 1, 1}, {0, 0, 0}), box({0.5f, 0.5f, 0.5f}, {0.8f, 1.4f, 0}), m));
    CHECK(m.count == 4);
    Real lo = 10, hi = -10;
    for (int i = 0; i < m.count; ++i) {
        lo = std::min(lo, m.points[i].point.x);
        hi = std::max(hi, m.points[i].point.x);
    }
    CHECK_NEAR(lo, 0.3, kEps);
    CHECK(hi >= 1.0f - 1e-4f && hi <= 1.0025f);  // clipped where the lower face ends (plus the 2 mm clip margin)
    CHECK_NEAR(patch_area_xz(m), 0.7, 1e-2);
}

TEST(boxes_are_separated_along_each_face_axis) {
    const Body a = box({1, 1, 1}, {0, 0, 0});
    Manifold m;
    for (int axis = 0; axis < 3; ++axis)
        for (Real sign : {-1.0f, 1.0f}) {
            Vec3 p;
            p[axis] = sign * 2.05f;
            CHECK(!collide(a, box({1, 1, 1}, p), m));
            p[axis] = sign * 1.95f;
            CHECK(collide(a, box({1, 1, 1}, p), m));
            Vec3 expected;
            expected[axis] = sign;
            CHECK_VEC(m.normal, expected, kEps);
            CHECK_NEAR(m.depth, 0.05, kEps);
        }
}

// Two cubes each turned 45 degrees, one about z and one about y, so that an edge of each points at the
// other and the two edges cross at right angles. For centres between 2.83 and 3.83 apart NO face normal of
// either cube separates them; only the direction perpendicular to both edges does. A separating axis test
// that checked just the six face normals would report a collision that is not there.
TEST(edge_against_edge_needs_the_cross_product_axes) {
    const Body a = box({1, 1, 1}, {0, 0, 0}, Quat::from_axis_angle({0, 0, 1}, kPi / 4));
    auto other = [](Real x) { return box({1, 1, 1}, {x, 0, 0}, Quat::from_axis_angle({0, 1, 0}, kPi / 4)); };
    Manifold m;
    CHECK(!collide(a, other(3.0f), m));  // edges 0.17 apart
    CHECK(!collide(a, other(2.9f), m));

    CHECK(collide(a, other(2.7f), m));   // edges overlap by 2 sqrt(2) - 2.7
    CHECK_VEC(m.normal, Vec3(1, 0, 0), kEps);
    CHECK_NEAR(m.depth, 2 * std::sqrt(2.0) - 2.7, kEps);
    CHECK(m.count == 1);
    CHECK_VEC(m.points[0].point, Vec3(1.35f, 0, 0), kEps);  // midway between the two edges, where they cross
}

TEST(a_corner_pressed_into_a_face_gives_one_point) {
    // Turn a cube so its body diagonal (1,1,1) points straight down, and lower that corner into a face.
    const Vec3 diagonal = Vec3(1, 1, 1).normalized(), down{0, -1, 0};
    const Quat q = Quat::from_axis_angle(cross(diagonal, down), std::acos(dot(diagonal, down)));
    const Real reach = 0.5f * std::sqrt(3.0f);  // centre to corner of a half-extent 0.5 cube
    Manifold m;
    CHECK(collide(box({2, 1, 2}, {0, 0, 0}), box({0.5f, 0.5f, 0.5f}, {0.3f, 1 + reach - 0.05f, -0.2f}, q), m));
    CHECK_VEC(m.normal, Vec3(0, 1, 0), kEps);
    CHECK_NEAR(m.depth, 0.05, kEps);
    CHECK(m.count == 1);
    CHECK_VEC(m.points[0].point, Vec3(0.3f, 0.975f, -0.2f), 2e-3);
}

TEST(eight_clipped_points_are_reduced_to_four_that_keep_most_of_the_patch) {
    // Equal cubes, the upper one turned 45 degrees about the vertical: the faces overlap in a regular
    // octagon (area 8 tan(22.5 deg) = 3.31). Four of its corners must be kept, and they should span a
    // good part of it; the best possible four (every other corner) cover 2.34.
    Manifold m;
    CHECK(collide(box({1, 1, 1}, {0, 0, 0}), box({1, 1, 1}, {0, 1.9f, 0}, Quat::from_axis_angle({0, 1, 0}, kPi / 4)), m));
    CHECK(m.count == 4);
    CHECK_NEAR(m.depth, 0.1, kEps);
    for (int i = 0; i < m.count; ++i) CHECK_NEAR(m.points[i].point.y, 0.95, kEps);
    CHECK(patch_area_xz(m) > 2.0);
    const auto ids = sorted_ids(m);
    CHECK(std::adjacent_find(ids.begin(), ids.end()) == ids.end());
}

TEST(contact_ids_stay_put_while_a_box_slides) {
    const Body floor = Body::fixed_box({5, 0.5f, 5}, {0, 0, 0});
    std::vector<std::uint32_t> first;
    for (int step = 0; step <= 30; ++step) {
        Manifold m;
        CHECK(collide(floor, box({0.5f, 0.5f, 0.5f}, {-1.5f + 0.1f * static_cast<Real>(step), 0.99f, 0.3f}, Quat::from_axis_angle({0, 1, 0}, 0.3f)), m));
        CHECK(m.count == 4);
        const auto ids = sorted_ids(m);
        if (step == 0) first = ids;
        CHECK(ids == first);
    }
    // Lying on a different face is a different set of features.
    Manifold other;
    CHECK(collide(floor, box({0.5f, 0.5f, 0.5f}, {0, 0.99f, 0.3f}, Quat::from_axis_angle({1, 0, 0}, kPi / 2)), other));
    CHECK(sorted_ids(other) != first);
}

// ---- properties over many random pairs ----------------------------------------------------------------

// The reported (normal, depth) must be a real way out: move B along the normal by a hair more than depth
// and the pair is disjoint; move it by clearly less and they still overlap.
TEST(random_pairs_report_a_translation_that_separates_them) {
    Lcg rng;
    int collisions = 0, box_box = 0;
    for (int trial = 0; trial < 4000; ++trial) {
        const Body a = random_body(rng), b = random_body(rng);
        Manifold m;
        if (!collide(a, b, m)) continue;
        ++collisions;
        box_box += a.shape.type == Shape::Type::Box && b.shape.type == Shape::Type::Box;

        CHECK_NEAR(m.normal.length(), 1, 1e-4);
        CHECK(m.count >= 1 && m.count <= Manifold::kMaxPoints);
        CHECK(m.depth >= 0);

        Body past = b;
        past.pos = b.pos + m.normal * (m.depth + 0.003f);
        Manifold unused;
        CHECK(!collide(a, past, unused));

        // The face-preference bias may report an axis up to 5% + 1 cm deeper than the true minimum, so
        // "clearly less" allows for that.
        const Real partial = 0.9f * m.depth - 0.03f;
        if (partial > 0.01f) {
            Body short_of = b;
            short_of.pos = b.pos + m.normal * partial;
            CHECK(collide(a, short_of, unused));
        }
    }
    CHECK(collisions > 1500);
    CHECK(box_box > 500);
}

// Every contact point lies in (or within the overlap depth of) both bodies, and its own depth is sane.
TEST(random_pairs_put_their_contact_points_in_the_overlap) {
    Lcg rng;
    rng.s = 77;
    int checked = 0;
    for (int trial = 0; trial < 3000; ++trial) {
        const Body a = random_body(rng), b = random_body(rng);
        Manifold m;
        if (!collide(a, b, m)) continue;
        for (int i = 0; i < m.count; ++i) {
            CHECK(contains_inflated(a, m.points[i].point, m.depth + 2e-3f));
            CHECK(contains_inflated(b, m.points[i].point, m.depth + 2e-3f));
            CHECK(m.points[i].depth >= -1e-4f && m.points[i].depth <= m.depth + 2e-3f);
            ++checked;
        }
    }
    CHECK(checked > 2000);
}

// If some point lies inside both bodies, the narrow phase must say they collide.
TEST(random_pairs_agree_with_point_sampling) {
    Lcg rng;
    rng.s = 5;
    int shared = 0;
    for (int trial = 0; trial < 600; ++trial) {
        const Body a = random_body(rng), b = random_body(rng);
        bool found = false;
        for (int k = 0; k < 1500 && !found; ++k) found = contains(b, point_inside(a, rng));
        Manifold m;
        if (found) {
            ++shared;
            CHECK(collide(a, b, m));
        }
    }
    CHECK(shared > 200);
}

// Turning and moving the whole scene must turn the normal with it and leave the depth alone.
TEST(collision_is_invariant_under_rigid_motion_in_3d) {
    Lcg rng;
    rng.s = 404;
    int checked = 0;
    for (int trial = 0; trial < 2500; ++trial) {
        const Body a = random_body(rng), b = random_body(rng);
        Manifold m1;
        if (!collide(a, b, m1)) continue;
        const Quat turn = rng.rotation();
        const Vec3 shift = rng.vec(6);
        Body a2 = a, b2 = b;
        a2.pos = rotate(turn, a.pos) + shift;
        b2.pos = rotate(turn, b.pos) + shift;
        a2.q = (turn * a.q).normalized();
        b2.q = (turn * b.q).normalized();
        Manifold m2;
        if (!collide(a2, b2, m2)) {
            CHECK(m1.depth < 2e-3f);  // only a razor-thin overlap may be lost to rounding
            continue;
        }
        // Two candidate axes within rounding of the bias threshold may be chosen differently after the
        // turn; then both answers are valid and close in depth. Otherwise the normals must match.
        if (dot(rotate(turn, m1.normal), m2.normal) > 0.999f) {
            CHECK_NEAR(m2.depth, m1.depth, 3e-3);
            ++checked;
        } else {
            CHECK_NEAR(m2.depth, m1.depth, 0.05f * m1.depth + 0.02f);
        }
    }
    CHECK(checked > 1000);
}

TEST(swapping_the_bodies_flips_the_normal) {
    Lcg rng;
    rng.s = 31;
    int checked = 0;
    for (int trial = 0; trial < 2000; ++trial) {
        const Body a = random_body(rng), b = random_body(rng);
        Manifold ab, ba;
        const bool hit = collide(a, b, ab);
        CHECK(collide(b, a, ba) == hit || ab.depth < 2e-3f);
        if (!hit || !collide(b, a, ba)) continue;
        // "A" keeps its own face unless another axis is clearly shallower, and after the swap a different
        // box is A. So the two calls may settle on different, nearly tied axes: their depths then differ by
        // at most the bias (5% + 1 cm). When they did pick the same axis, the depths match and the normals
        // must be exactly opposite.
        CHECK_NEAR(ab.depth, ba.depth, 0.06f * std::max(ab.depth, ba.depth) + 0.012f);
        if (std::fabs(ab.depth - ba.depth) < 1e-4f && ab.depth > 0.01f) {
            CHECK(dot(ab.normal, ba.normal) < -0.99f);
            ++checked;
        }
    }
    CHECK(checked > 700);
}

TEST(contains_tests_points_in_the_bodys_own_frame) {
    const Body tilted = box({2, 0.5f, 0.5f}, {1, 1, 1}, Quat::from_axis_angle({0, 0, 1}, kPi / 2));  // long axis now along y
    CHECK(contains(tilted, {1, 2.8f, 1}));
    CHECK(!contains(tilted, {2.8f, 1, 1}));
    CHECK(contains(ball(1, {0, 0, 0}), {0.5f, 0.5f, 0.5f}));
    CHECK(!contains(ball(1, {0, 0, 0}), {0.7f, 0.7f, 0.7f}));
}

// A box lying almost flat on another, at any heading and with a small random tilt, must still report a
// patch of at least three points and a near-vertical normal: enough to stop it rocking.
TEST(a_nearly_flat_contact_keeps_its_patch) {
    Lcg rng;
    rng.s = 99;
    const Body below = box({1, 1, 1}, {0, 0, 0});
    int patches = 0;
    for (int trial = 0; trial < 400; ++trial) {
        const Quat tilt = Quat::from_axis_angle(rng.vec(1), rng.range(0.0005f, 0.01f));
        const Quat turn = Quat::from_axis_angle({0, 1, 0}, rng.range(0, 6.28f));  // any heading
        Manifold m;
        CHECK(collide(below, box({0.6f, 0.4f, 0.5f}, {rng.range(-0.2f, 0.2f), 1.37f, rng.range(-0.2f, 0.2f)}, tilt * turn), m));
        CHECK(m.count >= 3);
        CHECK(m.normal.y > 0.999f);
        patches += m.count >= 3;
    }
    CHECK(patches == 400);
}

// The same face of a box, pressed against two different faces of another, makes two different contacts.
TEST(a_different_reference_face_gives_different_ids) {
    const Body a = box({1, 1, 1}, {0, 0, 0});
    Manifold on_top, at_side;
    CHECK(collide(a, box({0.5f, 0.5f, 0.5f}, {0, 1.45f, 0}), on_top));  // B's -y face on A's +y face
    // Turned a quarter turn the other way about z, B's own -y face looks toward -x: onto A's +x face.
    CHECK(collide(a, box({0.5f, 0.5f, 0.5f}, {1.45f, 0, 0}, Quat::from_axis_angle({0, 0, 1}, -kPi / 2)), at_side));
    CHECK_VEC(on_top.normal, Vec3(0, 1, 0), kEps);
    CHECK_VEC(at_side.normal, Vec3(1, 0, 0), kEps);
    CHECK(on_top.count == 4 && at_side.count == 4);
    for (std::uint32_t id : sorted_ids(on_top))
        for (std::uint32_t other : sorted_ids(at_side)) CHECK(id != other);
}

// Equal cubes stacked: the upper face of the lower cube and the lower face of the upper one are the same
// plane, so "A's face" and "B's face" overlap by the same amount up to rounding. Picking the strict
// minimum would make the reference face flip between the two boxes from one frame to the next as the
// stack jitters, and every flip changes all the contact ids (throwing away the warm start). The
// preference for A's face unless another is clearly shallower keeps it steady.
TEST(a_tied_face_pair_keeps_the_same_reference_under_jitter) {
    Lcg rng;
    rng.s = 1234;
    const Body below = box({0.5f, 0.5f, 0.5f}, {0, 0, 0});
    std::vector<std::uint32_t> first;
    for (int trial = 0; trial < 300; ++trial) {
        const Quat wobble = Quat::from_axis_angle(rng.vec(1), rng.range(0, 2e-4f));
        Manifold m;
        CHECK(collide(below, box({0.5f, 0.5f, 0.5f}, Vec3{0, 0.99f, 0} + rng.vec(1e-4f), wobble), m));
        CHECK(m.count == 4);
        const auto ids = sorted_ids(m);
        if (trial == 0) first = ids;
        CHECK(ids == first);
    }
}
