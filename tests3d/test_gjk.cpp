#include "test.hpp"
#include <phys3d/collision.hpp>
#include <phys3d/gjk.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace phys3d;

namespace {

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Vec3 in_box(Real r) { return {range(-r, r), range(-r, r), range(-r, r)}; }
    Vec3 direction() {
        for (;;) {
            const Vec3 v = in_box(1);
            if (v.length_sq() > 0.01f && v.length_sq() <= 1) return v.normalized();
        }
    }
    Quat rotation() { return Quat::from_axis_angle(direction(), range(-kPi, kPi)); }
};

// How far A and B overlap when both are flattened onto direction d: the width of A - B's shadow past the
// origin. Negative means d is a separating direction, with that much room. Uses nothing but the support
// functions, so it is an independent check on anything GJK and EPA claim.
Real overlap_along(const Convex& a, const Convex& b, Vec3 d) { return dot(d, a.support(d) - b.support(-d)); }

// Is p inside (or on) the shape, to within tol? A point is inside a convex shape exactly when no direction
// has it beyond the shape's farthest point.
bool inside(const Convex& shape, Vec3 p, Real tol, Lcg& rng, Vec3 also) {
    if (dot(p - shape.support(also), also) > tol || dot(p - shape.support(-also), -also) > tol) return false;
    for (int i = 0; i < 60; ++i) {
        const Vec3 d = rng.direction();
        if (dot(p - shape.support(d), d) > tol) return false;
    }
    return true;
}

// Checks a result against the definition of what it claims, using only support functions. Returns a
// description of the first thing wrong, or nullptr.
const char* audit(const Convex& a, const Convex& b, const ClosestResult& r, Real tol, Lcg& rng) {
    if (!std::isfinite(r.distance) || !std::isfinite(r.normal.x + r.normal.y + r.normal.z)) return "not finite";
    if (std::fabs(r.normal.length() - 1) > 1e-3f) return "normal is not a unit vector";
    if (!inside(a, r.point_a, tol, rng, r.normal)) return "point_a is not on a";
    if (!inside(b, r.point_b, tol, rng, r.normal)) return "point_b is not on b";
    // The two points are `distance` apart along the normal (behind one another when overlapping).
    if (distance(r.point_a + r.normal * r.distance, r.point_b) > tol) return "points do not match distance and normal";
    // Flattened onto the normal, the shapes are apart by exactly `distance` (overlap by exactly the depth):
    // moving b that far along the normal separates them, so the true answer is no bigger...
    if (std::fabs(overlap_along(a, b, r.normal) + r.distance) > tol) return "normal does not give the claimed distance";
    if (r.distance < 0) {
        // ...and for an overlapping pair no other direction may get out sooner. Directions close to the
        // normal are the ones that could.
        for (int i = 0; i < 80; ++i) {
            const Vec3 d = i < 40 ? rng.direction() : (r.normal + rng.direction() * rng.range(0.001f, 0.15f)).normalized();
            if (overlap_along(a, b, d) < -r.distance - tol) return "a shorter way out exists";
        }
    }
    return nullptr;
}

// The exact penetration depth of two boxes by the separating axis theorem: the least overlap over the 15
// candidate axes. (collide() has this too, but biased toward face axes; this is the plain minimum.)
Real box_depth_by_sat(const Convex& a, const Convex& b) {
    const Mat3 ra = to_mat3(a.q), rb = to_mat3(b.q);
    const Vec3 ea[3] = {ra.cx, ra.cy, ra.cz}, eb[3] = {rb.cx, rb.cy, rb.cz};
    std::vector<Vec3> axes(ea, ea + 3);
    axes.insert(axes.end(), eb, eb + 3);
    for (const Vec3& x : ea)
        for (const Vec3& y : eb)
            if (cross(x, y).length() > 1e-3f) axes.push_back(cross(x, y).normalized());
    Real least = 1e30f;
    for (const Vec3& d : axes) least = std::min({least, overlap_along(a, b, d), overlap_along(a, b, -d)});
    return least;
}

struct RandomShape {
    std::vector<Vec3> points;  // storage for a hull
    Convex convex;
};
// kind: 0 sphere, 1 box, 2 capsule, 3 hull, 4 rounded hull.
void make_shape(RandomShape& out, int kind, Lcg& rng, Vec3 pos, Real scale) {
    const Quat q = rng.rotation();
    switch (kind) {
        case 0: out.convex = Convex::sphere(rng.range(0.1f, 1) * scale, pos); break;
        case 1: out.convex = Convex::box(Vec3{rng.range(0.1f, 1), rng.range(0.1f, 1), rng.range(0.1f, 1)} * scale, pos, q); break;
        case 2: out.convex = Convex::capsule(rng.range(0.1f, 1) * scale, rng.range(0.05f, 0.5f) * scale, pos, q); break;
        default: {
            const int n = 4 + static_cast<int>(rng.next() * 9);
            out.points.clear();
            for (int i = 0; i < n; ++i) out.points.push_back(rng.in_box(1) * scale);
            out.convex = Convex::hull(out.points.data(), n, pos, q, kind == 4 ? 0.15f * scale : Real(0));
        }
    }
}

}  // namespace

// ---- the two routines underneath -----------------------------------------------------------------

namespace {

// The distance from the origin to triangle abc, found the slow way: the foot of the perpendicular if it
// lands inside, otherwise the nearest point of the three edges, each by its own clamped projection.
Real triangle_distance_slowly(Vec3 a, Vec3 b, Vec3 c) {
    auto edge = [](Vec3 p, Vec3 q) {
        const Vec3 pq = q - p;
        const Real t = std::clamp(dot(-p, pq) / std::max(pq.length_sq(), 1e-30f), 0.0f, 1.0f);
        return (p + pq * t).length();
    };
    Real best = std::min({edge(a, b), edge(b, c), edge(c, a)});
    const Vec3 n = cross(b - a, c - a);
    if (n.length() > 1e-6f) {
        const Vec3 unit = n.normalized(), foot = unit * dot(unit, a);
        const bool inside = dot(cross(b - a, foot - a), n) >= 0 && dot(cross(c - b, foot - b), n) >= 0 && dot(cross(a - c, foot - c), n) >= 0;
        if (inside) best = std::min(best, foot.length());
    }
    return best;
}

}  // namespace

TEST(closest_point_on_a_segment_is_right_at_both_ends_and_between) {
    Real w[2];
    detail::closest_on_segment({1, 2, 0}, {3, 2, 0}, w);  // the origin is behind a
    CHECK(w[0] == 1 && w[1] == 0);
    detail::closest_on_segment({-3, 2, 0}, {-1, 2, 0}, w);  // beyond b
    CHECK(w[0] == 0 && w[1] == 1);
    detail::closest_on_segment({-1, 2, 0}, {3, 2, 0}, w);  // a quarter of the way along
    CHECK_NEAR(w[0], 0.75, 1e-6);
    CHECK_NEAR(w[1], 0.25, 1e-6);
    detail::closest_on_segment({1, 1, 1}, {1, 1, 1}, w);  // no length at all
    CHECK_NEAR(w[0] + w[1], 1, 1e-6);
}

TEST(closest_point_on_a_triangle_is_right_in_all_seven_regions) {
    Lcg rng;
    int corners = 0, edges = 0, inside = 0;
    for (int trial = 0; trial < 20000; ++trial) {
        const Vec3 a = rng.in_box(1), b = rng.in_box(1), c = rng.in_box(1);
        if (cross(b - a, c - a).length() < 0.05f) continue;
        Real w[3];
        detail::closest_on_triangle(a, b, c, w);
        CHECK(w[0] >= 0 && w[1] >= 0 && w[2] >= 0);
        CHECK_NEAR(w[0] + w[1] + w[2], 1, 1e-5);
        CHECK_NEAR((a * w[0] + b * w[1] + c * w[2]).length(), triangle_distance_slowly(a, b, c), 1e-4);
        const int zeros = (w[0] == 0) + (w[1] == 0) + (w[2] == 0);
        (zeros == 2 ? corners : zeros == 1 ? edges : inside) += 1;
        // Which corner or edge it lands on must not depend on the order the corners are given in.
        Real v[3];
        detail::closest_on_triangle(c, a, b, v);
        CHECK(distance(a * w[0] + b * w[1] + c * w[2], c * v[0] + a * v[1] + b * v[2]) < 1e-4f);
    }
    CHECK(corners > 400 && edges > 2000 && inside > 2000);  // every kind of region was exercised
}

TEST(a_triangle_with_no_area_still_gets_an_answer) {
    // Three points in a line, the origin off to one side of the middle stretch.
    Real w[3];
    detail::closest_on_triangle({-2, 1, 0}, {4, 1, 0}, {1, 1, 0}, w);
    CHECK(std::isfinite(w[0] + w[1] + w[2]));
    CHECK_NEAR(w[0] + w[1] + w[2], 1, 1e-5);
    CHECK(distance(Vec3{-2, 1, 0} * w[0] + Vec3{4, 1, 0} * w[1] + Vec3{1, 1, 0} * w[2], Vec3{0, 1, 0}) < 1e-5f);
    // And all three in one place.
    detail::closest_on_triangle({1, 2, 3}, {1, 2, 3}, {1, 2, 3}, w);
    CHECK_NEAR(w[0] + w[1] + w[2], 1, 1e-5);
}

// ---- known answers: apart ------------------------------------------------------------------------

TEST(gjk_gives_the_distance_between_two_spheres) {
    const ClosestResult r = closest(Convex::sphere(0.5f, {1, 2, 3}), Convex::sphere(0.25f, {1, 2, 5}));
    CHECK_NEAR(r.distance, 2 - 0.75, 1e-5);
    CHECK(distance(r.normal, Vec3{0, 0, 1}) < 1e-5f);
    CHECK(distance(r.point_a, Vec3{1, 2, 3.5f}) < 1e-5f);
    CHECK(distance(r.point_b, Vec3{1, 2, 4.75f}) < 1e-5f);
    CHECK(r.epa_iterations == 0);
}

TEST(gjk_gives_the_gap_between_two_boxes_face_to_face) {
    const ClosestResult r = closest(Convex::box({1, 1, 1}, {0, 0, 0}), Convex::box({0.5f, 2, 0.5f}, {2, 0.3f, 0.2f}));
    CHECK_NEAR(r.distance, 0.5, 1e-5);
    CHECK(distance(r.normal, Vec3{1, 0, 0}) < 1e-5f);
    CHECK_NEAR(r.point_a.x, 1, 1e-5);
    CHECK_NEAR(r.point_b.x, 1.5, 1e-5);
    CHECK(distance(r.point_a + Vec3{0.5f, 0, 0}, r.point_b) < 1e-5f);  // directly opposite each other
}

TEST(gjk_gives_the_distance_from_a_box_corner_to_a_face) {
    // A cube balanced on one corner (its long diagonal vertical) above a slab.
    const Vec3 diagonal = Vec3{1, 1, 1}.normalized();
    const Vec3 axis = cross(diagonal, Vec3{0, -1, 0});
    const Quat q = Quat::from_axis_angle(axis, std::acos(dot(diagonal, Vec3{0, -1, 0})));
    const Real corner_drop = std::sqrt(3.0f) * 0.5f;  // centre to corner of a unit cube
    const ClosestResult r = closest(Convex::box({0.5f, 0.5f, 0.5f}, {0.3f, 2, -0.2f}, q), Convex::box({3, 0.5f, 3}, {0, 0, 0}));
    CHECK_NEAR(r.distance, 2 - corner_drop - 0.5, 1e-4);
    CHECK(distance(r.normal, Vec3{0, -1, 0}) < 1e-4f);
    CHECK(distance(r.point_a, Vec3{0.3f, 2 - corner_drop, -0.2f}) < 1e-4f);
    CHECK(distance(r.point_b, Vec3{0.3f, 0.5f, -0.2f}) < 1e-4f);
}

TEST(gjk_gives_the_distance_between_two_crossed_edges) {
    // Two long thin bars, one along x and one along z a little higher, each turned 45 degrees about its own
    // length so that an edge, not a face, is what faces the other bar.
    const Real t = 0.1f, edge = t * std::sqrt(2.0f);
    const ClosestResult r = closest(Convex::box({2, t, t}, {0, 0, 0}, Quat::from_axis_angle({1, 0, 0}, kPi / 4)),
                                    Convex::box({t, t, 2}, {0.4f, 1, -0.3f}, Quat::from_axis_angle({0, 0, 1}, kPi / 4)));
    CHECK_NEAR(r.distance, 1 - 2 * edge, 1e-4);
    CHECK(distance(r.normal, Vec3{0, 1, 0}) < 1e-3f);
    CHECK(distance(r.point_a, Vec3{0.4f, edge, 0}) < 1e-3f);
    CHECK(distance(r.point_b, Vec3{0.4f, 1 - edge, 0}) < 1e-3f);
}

TEST(gjk_gives_the_distance_between_capsules) {
    // Crossed at right angles, 3 apart: the gap is 3 less the two radii.
    const Convex upright = Convex::capsule(1, 0.3f, {0, 0, 0});
    const Convex lying = Convex::capsule(1, 0.5f, {3, 0.2f, 0}, Quat::from_axis_angle({1, 0, 0}, kPi / 2));
    const ClosestResult r = closest(upright, lying);
    CHECK_NEAR(r.distance, 3 - 0.8, 1e-5);
    CHECK(distance(r.normal, Vec3{1, 0, 0}) < 1e-5f);
    // End to end along y: the caps face each other.
    const ClosestResult ends = closest(upright, Convex::capsule(0.5f, 0.2f, {0, 4, 0}));
    CHECK_NEAR(ends.distance, 4 - 1 - 0.5 - 0.3 - 0.2, 1e-5);
}

TEST(gjk_handles_a_hull_exactly_like_the_box_it_is) {
    const Vec3 h{0.4f, 0.7f, 0.2f};
    Vec3 corners[8];
    for (int i = 0; i < 8; ++i) corners[i] = {(i & 1) ? h.x : -h.x, (i & 2) ? h.y : -h.y, (i & 4) ? h.z : -h.z};
    Lcg rng;
    for (int trial = 0; trial < 300; ++trial) {
        const Vec3 pa = rng.in_box(1), pb = rng.in_box(1.5f);
        const Quat qa = rng.rotation(), qb = rng.rotation();
        const Convex other = Convex::box({0.5f, 0.3f, 0.6f}, pb, qb);
        const ClosestResult as_box = closest(Convex::box(h, pa, qa), other);
        const ClosestResult as_hull = closest(Convex::hull(corners, 8, pa, qa), other);
        CHECK_NEAR(as_hull.distance, as_box.distance, 1e-3);
    }
}

// ---- known answers: overlapping ------------------------------------------------------------------

TEST(overlapping_spheres_need_no_epa_and_are_exact) {
    const ClosestResult r = closest(Convex::sphere(1, {0, 0, 0}), Convex::sphere(0.5f, {0, 1.2f, 0}));
    CHECK_NEAR(r.distance, -0.3, 1e-6);
    CHECK(distance(r.normal, Vec3{0, 1, 0}) < 1e-6f);
    CHECK(distance(r.point_a, Vec3{0, 1, 0}) < 1e-6f);     // a's surface, inside b
    CHECK(distance(r.point_b, Vec3{0, 0.7f, 0}) < 1e-6f);  // b's surface, inside a
    CHECK(r.epa_iterations == 0);
}

TEST(epa_gives_the_depth_of_a_small_box_inside_a_big_one) {
    // The shortest way out is through the nearest face: (big half + small half - offset) on that axis.
    const ClosestResult r = closest(Convex::box({2, 3, 4}, {0, 0, 0}), Convex::box({0.2f, 0.2f, 0.2f}, {1.5f, 0.5f, -1}));
    CHECK_NEAR(r.distance, -(2 + 0.2 - 1.5), 1e-4);
    CHECK(distance(r.normal, Vec3{1, 0, 0}) < 1e-4f);
    CHECK(r.epa_iterations > 0);
}

TEST(epa_gives_the_depth_of_a_point_inside_a_box) {
    const ClosestResult r = closest(Convex::box({1, 2, 3}, {5, 5, 5}), Convex::sphere(0, {5.2f, 6.7f, 4}));
    CHECK_NEAR(r.distance, -0.3, 1e-4);  // 2 - 1.7 through the top
    CHECK(distance(r.normal, Vec3{0, 1, 0}) < 1e-4f);
    CHECK(distance(r.point_a, Vec3{5.2f, 7, 4}) < 1e-3f);
}

TEST(two_identical_boxes_in_the_same_place_come_apart_along_the_thinnest_axis) {
    const Quat q = Quat::from_axis_angle({1, 2, 3}, 0.8f);
    const Convex box = Convex::box({0.6f, 0.2f, 0.9f}, {1, 1, 1}, q);
    const ClosestResult r = closest(box, box);
    CHECK_NEAR(r.distance, -0.4, 1e-4);
    CHECK_NEAR(std::fabs(dot(r.normal, rotate(q, {0, 1, 0}))), 1, 1e-4);
}

TEST(concentric_spheres_come_apart_by_the_sum_of_their_radii) {
    const ClosestResult r = closest(Convex::sphere(0.7f, {2, 2, 2}), Convex::sphere(0.4f, {2, 2, 2}));
    CHECK_NEAR(r.distance, -1.1, 1e-6);
    CHECK_NEAR(r.normal.length(), 1, 1e-6);  // any direction will do, but it must be one
}

TEST(capsules_whose_axes_cross_come_apart_perpendicular_to_both) {
    // The cores are two segments meeting at a point: their Minkowski difference is a flat parallelogram, so
    // there is no tetrahedron for EPA to grow. The way out is along the parallelogram's normal.
    const Convex a = Convex::capsule(1, 0.3f, {0, 0, 0});                                              // along y
    const Convex b = Convex::capsule(1, 0.2f, {0.1f, 0.2f, 0}, Quat::from_axis_angle({0, 0, 1}, 1.0f));  // in the xy plane
    const ClosestResult r = closest(a, b);
    CHECK_NEAR(r.distance, -0.5, 1e-4);
    CHECK_NEAR(std::fabs(r.normal.z), 1, 1e-4);
    CHECK(distance(r.point_b + r.normal * 0.5f, r.point_a) < 1e-4f);
}

TEST(capsules_that_pass_a_hair_apart_still_get_an_accurate_normal) {
    // Found by a million-pair stress run. When two axes miss each other by a millimetre or less, the
    // closest points are a tiny difference of large numbers, and a normal taken from them was a third of
    // a degree out: enough to misjudge the overlap by 3.5 mm. The common perpendicular of the two axes is
    // the only possible normal.
    Lcg rng;
    Real worst_angle = 0, worst_depth = 0;
    for (int trial = 0; trial < 300; ++trial) {
        const Quat qa = rng.rotation(), qb = rng.rotation();
        const Vec3 axis_a = rotate(qa, {0, 1, 0}), axis_b = rotate(qb, {0, 1, 0});
        const Vec3 perpendicular = cross(axis_a, axis_b);
        if (perpendicular.length() < 0.3f) continue;  // nearly parallel axes are a different case
        const Vec3 n = perpendicular.normalized();
        // Each axis passes through a point near its own middle; the two points are `miss` apart along n.
        const Vec3 crossing = rng.in_box(0.5f);
        const Real miss = rng.range(5e-4f, 2e-3f);  // (anything under a ten-thousandth of the size counts as touching)
        const Convex a = Convex::capsule(1, 0.2f, crossing + axis_a * rng.range(-0.6f, 0.6f), qa);
        const Convex b = Convex::capsule(1, 0.3f, crossing + n * miss + axis_b * rng.range(-0.6f, 0.6f), qb);
        const ClosestResult r = closest(a, b);
        worst_angle = std::max(worst_angle, cross(r.normal, n).length());
        worst_depth = std::max(worst_depth, std::fabs(-r.distance - (0.5f - miss)));
    }
    CHECK(worst_angle < 2e-4f);
    CHECK(worst_depth < 2e-5f);
}

TEST(a_point_exactly_on_a_corner_edge_or_face_of_a_box_is_handled) {
    // GJK's very first sample of A - B is then the origin itself, and EPA has to build its starting
    // tetrahedron from a single point.
    Lcg rng;
    const Convex box = Convex::box({1, 1, 1}, {0, 0, 0});
    // (Off-centre points too: for those the direction between the two centres is not a way out.)
    for (const Vec3& on : {Vec3{1, 1, 1}, Vec3{1, 1, 0}, Vec3{1, 0, 0}, Vec3{-1, 1, -1}, Vec3{0, -1, 1}, Vec3{1, 0.5f, 0.3f}, Vec3{0.2f, -1, 0.7f}, Vec3{1, 1, 0.4f}, Vec3{-0.6f, 1, -1}}) {
        const Convex ball = Convex::sphere(0.25f, on);
        const ClosestResult r = closest(box, ball);
        CHECK(r.converged);
        CHECK_NEAR(r.distance, -0.25, 1e-4);
        const char* wrong = audit(box, ball, r, 1e-3f, rng);
        if (wrong) std::printf("    at (%g, %g, %g): %s\n", on.x, on.y, on.z, wrong);
        CHECK(wrong == nullptr);
        const ClosestResult swapped = closest(ball, box);
        CHECK_NEAR(swapped.distance, -0.25, 1e-4);
        CHECK(audit(ball, box, swapped, 1e-3f, rng) == nullptr);
    }
}

TEST(touching_shapes_report_zero) {
    const ClosestResult faces = closest(Convex::box({1, 1, 1}, {0, 0, 0}), Convex::box({1, 1, 1}, {2, 0.5f, 0}));
    CHECK_NEAR(faces.distance, 0, 1e-4);
    CHECK_NEAR(faces.normal.x, 1, 1e-3);  // from the first box toward the second, not the other way
    const ClosestResult balls = closest(Convex::sphere(1, {0, 0, 0}), Convex::sphere(2, {0, 3, 0}));
    CHECK_NEAR(balls.distance, 0, 1e-6);
}

// ---- against the code that knows the shapes -------------------------------------------------------

TEST(epa_agrees_with_the_separating_axis_theorem_for_boxes) {
    Lcg rng;
    int overlapping = 0;
    Real worst = 0;
    for (int trial = 0; trial < 3000; ++trial) {
        const Convex a = Convex::box({rng.range(0.1f, 1), rng.range(0.1f, 1), rng.range(0.1f, 1)}, rng.in_box(1), rng.rotation());
        const Convex b = Convex::box({rng.range(0.1f, 1), rng.range(0.1f, 1), rng.range(0.1f, 1)}, rng.in_box(1), rng.rotation());
        const Real exact = box_depth_by_sat(a, b);
        const ClosestResult r = closest(a, b);
        CHECK((exact > 0) == (r.distance < 0) || std::fabs(exact) < 1e-4f);
        if (exact <= 0) continue;
        ++overlapping;
        worst = std::max(worst, std::fabs(-r.distance - exact));
    }
    CHECK(overlapping > 1500);
    CHECK(worst < 1e-3f);
}

TEST(gjk_and_epa_agree_with_collide_for_a_sphere_and_a_box) {
    Lcg rng;
    int hits = 0;
    for (int trial = 0; trial < 2000; ++trial) {
        const Body box = Body::solid_box({rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)}, 1, rng.in_box(1), rng.rotation());
        const Body ball = Body::solid_sphere(rng.range(0.1f, 0.8f), 1, rng.in_box(1.5f));
        Manifold m;
        const bool hit = collide(box, ball, m);
        const ClosestResult r = closest(Convex::of(box), Convex::of(ball));
        CHECK(hit == (r.distance < 0) || std::fabs(r.distance) < 1e-4f);
        if (!hit || r.distance >= 0) continue;
        ++hits;
        CHECK_NEAR(-r.distance, m.depth, 1e-3);
        // The normals agree unless two ways out tie, which only happens for a centre deep inside the box.
        if (r.epa_iterations == 0) CHECK(dot(r.normal, m.normal) > 0.999f);
    }
    CHECK(hits > 500);
}

// ---- every pairing, checked against the definition -----------------------------------------------

TEST(every_pairing_of_shapes_passes_an_audit_by_support_functions) {
    Lcg rng;
    int apart = 0, shallow = 0, deep = 0, failures = 0;
    for (int trial = 0; trial < 6000; ++trial) {
        RandomShape a, b;
        make_shape(a, trial % 5, rng, rng.in_box(0.5f), 1);
        make_shape(b, (trial / 5) % 5, rng, rng.in_box(trial % 2 ? 2.2f : 0.8f), 1);  // half far, half close
        const ClosestResult r = closest(a.convex, b.convex);
        const char* wrong = !r.converged ? "did not converge" : audit(a.convex, b.convex, r, 2e-3f, rng);
        if (wrong && ++failures <= 5)
            std::printf("    trial %d (kinds %d, %d): %s; distance %.5f after %d + %d iterations\n", trial, trial % 5, (trial / 5) % 5,
                        wrong, r.distance, r.gjk_iterations, r.epa_iterations);
        (r.distance >= 0 ? apart : r.epa_iterations == 0 ? shallow : deep) += 1;
    }
    CHECK(failures == 0);
    CHECK(apart > 1000 && shallow > 200 && deep > 1000);  // all three paths were exercised
}

TEST(shapes_moved_exactly_into_contact_pass_the_same_audit) {
    // Random pairs are never exactly touching, and exactly touching is where a sign can go wrong: the
    // origin sits on the surface of A - B, on neither side of it. So take each pair, slide b along the
    // normal until the two just touch, and ask again.
    Lcg rng;
    int failures = 0;
    for (int trial = 0; trial < 4000; ++trial) {
        RandomShape a, b;
        make_shape(a, trial % 5, rng, rng.in_box(0.5f), 1);
        make_shape(b, (trial / 5) % 5, rng, rng.in_box(trial % 2 ? 2.2f : 0.8f), 1);
        const ClosestResult before = closest(a.convex, b.convex);
        b.convex.pos = b.convex.pos - before.normal * before.distance;
        const ClosestResult r = closest(a.convex, b.convex);
        const char* wrong = !r.converged ? "did not converge" : std::fabs(r.distance) > 1e-3f ? "not touching" : audit(a.convex, b.convex, r, 2e-3f, rng);
        if (wrong && ++failures <= 5)
            std::printf("    trial %d (kinds %d, %d): %s; distance %.5f after %d + %d iterations\n", trial, trial % 5, (trial / 5) % 5,
                        wrong, r.distance, r.gjk_iterations, r.epa_iterations);
    }
    CHECK(failures == 0);
}

TEST(a_rod_touching_a_box_corner_end_on_gets_a_real_way_out) {
    // The end of a capsule's axis sits exactly on a corner of a box, and the capsule leans away so that
    // the direction between the two centres is NOT a way to separate them. The cores meet at one point, so
    // every sample of A - B that matters is the origin, and EPA has to grow its tetrahedron from that.
    Lcg rng;
    const Convex box = Convex::box({1, 1, 1}, {0, 0, 0});
    for (const Vec3& lean : {Vec3{1, -1, 0.2f}, Vec3{0.1f, 1, -1}, Vec3{-1, 0.3f, 1}, Vec3{1, 1, 1}}) {
        const Vec3 along = lean.normalized(), corner{1, 1, 1};
        const Quat q = Quat::from_axis_angle(cross({0, 1, 0}, along), std::acos(along.y));
        const Convex rod = Convex::capsule(5, 0.1f, corner + rotate(q, {0, 5, 0}), q);
        CHECK(distance(rod.core_support(-along), corner) < 1e-5f);  // the axis really does end on the corner
        for (int order = 0; order < 2; ++order) {
            const Convex& a = order ? rod : box;
            const Convex& b = order ? box : rod;
            const ClosestResult r = closest(a, b);
            CHECK(r.converged);
            CHECK_NEAR(r.distance, -0.1, 1e-4);
            const char* wrong = audit(a, b, r, 1e-3f, rng);
            if (wrong) std::printf("    lean (%g, %g, %g), order %d: %s\n", lean.x, lean.y, lean.z, order, wrong);
            CHECK(wrong == nullptr);
        }
    }
}

TEST(shapes_a_hair_from_contact_either_way_pass_the_same_audit) {
    // Slid to within two ten-thousandths of their size of touching, apart or overlapping. This is the band
    // where the closest points are too close together for their difference to give a direction.
    int failures = 0;
    for (int trial = 0; trial < 6000; ++trial) {
        Lcg rng;
        rng.s = static_cast<std::uint32_t>(trial) * 2654435761u + 17;
        const Real scale = trial % 3 == 0 ? 50.0f : 1.0f;
        RandomShape a, b;
        make_shape(a, trial % 5, rng, rng.in_box(0.5f) * scale, scale);
        make_shape(b, (trial / 5) % 5, rng, rng.in_box(trial % 2 ? 2.2f : 0.8f) * scale, scale);
        const ClosestResult before = closest(a.convex, b.convex);
        b.convex.pos = b.convex.pos - before.normal * (before.distance + rng.range(-2e-4f, 2e-4f) * scale);
        const ClosestResult r = closest(a.convex, b.convex);
        const char* wrong = !r.converged ? "did not converge" : audit(a.convex, b.convex, r, 2e-3f * scale, rng);
        if (wrong && ++failures <= 5)
            std::printf("    trial %d (kinds %d, %d): %s; distance %.6f after %d + %d iterations\n", trial, trial % 5, (trial / 5) % 5,
                        wrong, r.distance / scale, r.gjk_iterations, r.epa_iterations);
    }
    CHECK(failures == 0);
}

TEST(deepest_points_that_cannot_be_trusted_are_reported_as_not_converged) {
    // The one pair in two million from the stress run where the second GJK pass (which finds the deepest
    // points) stops early on rounding. The depth and normal are right; the points are about 5 mm out, and
    // the result has to say so.
    const Vec3 points[12] = {
        {0x1.e32a5p-1f, -0x1.916988p-1f, -0x1.8baa8p-1f}, {0x1.610cbp-2f, -0x1.de308p-1f, 0x1.1fb048p-1f},
        {-0x1.bfb008p-2f, 0x1.3f1bfp-1f, 0x1.e69c5p-2f},  {0x1.7a168p-2f, -0x1.a2f8a8p-1f, 0x1.17fbep-3f},
        {-0x1.a0d6a4p-1f, 0x1.c2bcap-3f, -0x1.171b7p-3f}, {0x1.a62f7p-2f, 0x1.e834ap-3f, 0x1.65b8d8p-1f},
        {0x1.989f8p-4f, 0x1.b95c6p-4f, -0x1.a10b2p-4f},   {-0x1.a8d5ep-1f, -0x1.60e16p-2f, -0x1.f2d348p-1f},
        {0x1.8ecbc8p-1f, 0x1.23b9c8p-1f, 0x1.0ba7ap-2f},  {0x1.8e8a9p-1f, 0x1.aeb78p-2f, 0x1.b51d4p-4f},
        {0x1.b32cp-3f, -0x1.7b7e98p-1f, 0x1.6afd18p-2f},  {-0x1.6a8774p-1f, -0x1.b0509p-1f, 0x1.266a8cp-1f}};
    const Convex a = Convex::hull(points, 12, {-0x1.1c3f8p-4f, 0x1.871bdp-4f, 0x1.2f4bf8p-2f},
                                  {-0x1.4728dep-1f, 0x1.458dap-3f, -0x1.f3eeb6p-2f, 0x1.25436ap-1f});
    const Convex b = Convex::box({0x1.d61c7ep-1f, 0x1.b72f0ep-1f, 0x1.18b66ap-3f}, {-0x1.4de896p-1f, -0x1.4a7496p-1f, -0x1.ef484ep-4f},
                                 {0x1.5c7b58p-3f, -0x1.43402ep-4f, -0x1.856006p-3f, 0x1.ed66eap-1f});
    const ClosestResult r = closest(a, b);
    CHECK_NEAR(overlap_along(a, b, r.normal), -r.distance, 1e-4);  // depth and normal agree with each other
    CHECK_NEAR(r.distance, -0.236254, 1e-4);
    Lcg rng;
    CHECK(!r.converged || audit(a, b, r, 2e-3f, rng) == nullptr);  // either the points are right, or it owns up
}

TEST(the_answer_does_not_depend_on_the_order_of_the_shapes) {
    Lcg rng;
    for (int trial = 0; trial < 1500; ++trial) {
        RandomShape a, b;
        make_shape(a, trial % 5, rng, rng.in_box(0.5f), 1);
        make_shape(b, (trial / 5) % 5, rng, rng.in_box(2), 1);
        const ClosestResult ab = closest(a.convex, b.convex), ba = closest(b.convex, a.convex);
        CHECK_NEAR(ab.distance, ba.distance, 2e-3);
    }
}

TEST(the_answer_does_not_depend_on_where_in_the_world_the_pair_is_or_how_big) {
    Lcg rng;
    Real worst_far = 0, worst_small = 0, worst_big = 0;
    for (int trial = 0; trial < 1500; ++trial) {
        const std::uint32_t seed = rng.s;
        auto distance_with = [&](Vec3 offset, Real scale) {
            rng.s = seed;
            RandomShape a, b;
            make_shape(a, trial % 5, rng, offset + rng.in_box(0.5f) * scale, scale);
            make_shape(b, (trial / 5) % 5, rng, offset + rng.in_box(2) * scale, scale);
            return closest(a.convex, b.convex).distance / scale;
        };
        const Real here = distance_with({0, 0, 0}, 1);
        // Positions a kilometre out are only stored to about 0.1 mm in single precision, so that is the
        // most that can be asked; the point is that the error is not far larger.
        worst_far = std::max(worst_far, std::fabs(distance_with({1000, -700, 400}, 1) - here));
        worst_small = std::max(worst_small, std::fabs(distance_with({0, 0, 0}, 0.01f) - here));
        worst_big = std::max(worst_big, std::fabs(distance_with({0, 0, 0}, 100) - here));
    }
    CHECK(worst_far < 3e-3f);
    CHECK(worst_small < 2e-3f);
    CHECK(worst_big < 2e-3f);
}

TEST(gjk_and_epa_finish_in_a_handful_of_steps) {
    Lcg rng;
    int worst_gjk = 0, worst_epa = 0;
    long total_gjk = 0;
    const int trials = 4000;
    for (int trial = 0; trial < trials; ++trial) {
        RandomShape a, b;
        make_shape(a, trial % 5, rng, rng.in_box(0.5f), 1);
        make_shape(b, (trial / 5) % 5, rng, rng.in_box(trial % 2 ? 2.2f : 0.8f), 1);
        const ClosestResult r = closest(a.convex, b.convex);
        worst_gjk = std::max(worst_gjk, r.gjk_iterations);
        worst_epa = std::max(worst_epa, r.epa_iterations);
        total_gjk += r.gjk_iterations;
    }
    CHECK(worst_gjk <= 30);
    CHECK(worst_epa <= 60);
    CHECK(total_gjk < 8L * trials);  // a handful on average
}
