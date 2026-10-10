#include "test.hpp"
#include <phys3d/collision.hpp>

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

Hull built(const std::vector<Vec3>& points) {
    Hull h;
    CHECK(build_hull(points.data(), static_cast<int>(points.size()), h));
    return h;
}
Hull box_hull(Vec3 half) {
    std::vector<Vec3> p;
    for (int i = 0; i < 8; ++i) p.push_back({(i & 1) ? half.x : -half.x, (i & 2) ? half.y : -half.y, (i & 4) ? half.z : -half.z});
    return built(p);
}
// A prism with a regular polygon of `sides` sides (circumradius r) for its two ends, along y.
Hull prism(int sides, Real r, Real half_height) {
    std::vector<Vec3> p;
    for (int i = 0; i < sides; ++i) {
        const Real a = 2 * kPi * static_cast<Real>(i) / static_cast<Real>(sides);
        p.push_back({r * std::cos(a), -half_height, r * std::sin(a)});
        p.push_back({r * std::cos(a), half_height, r * std::sin(a)});
    }
    return built(p);
}
Hull random_hull(Lcg& rng) {
    std::vector<Vec3> p;
    const int n = 6 + static_cast<int>(rng.next() * 14);
    for (int i = 0; i < n; ++i) p.push_back(rng.in_box(0.6f));
    return built(p);
}

int find_point(const Manifold& m, Vec3 p, Real tol) {
    for (int i = 0; i < m.count; ++i)
        if (distance(m.points[i].point, p) < tol) return i;
    return -1;
}
bool same_points(const Manifold& x, const Manifold& y, Real tol) {
    if (x.count != y.count) return false;
    for (int i = 0; i < x.count; ++i) {
        const int j = find_point(y, x.points[i].point, tol);
        if (j < 0 || std::fabs(y.points[j].depth - x.points[i].depth) > tol) return false;
    }
    return true;
}
bool inside(const Convex& shape, Vec3 p, Real tol, Lcg& rng, Vec3 also) {
    if (dot(p - shape.support(also), also) > tol || dot(p - shape.support(-also), -also) > tol) return false;
    for (int i = 0; i < 40; ++i) {
        const Vec3 d = rng.direction();
        if (dot(p - shape.support(d), d) > tol) return false;
    }
    return true;
}

const Convex kSlab = Convex::box({5, 0.5f, 5}, {0, -0.5f, 0});  // top face at y = 0

}  // namespace

// ---- a hull that is a box must behave as the box --------------------------------------------------

TEST(a_hull_built_from_a_boxs_corners_collides_exactly_as_the_box) {
    Lcg rng;
    int hits = 0, patches = 0;
    for (int trial = 0; trial < 4000; ++trial) {
        const Vec3 half_a{rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)}, half_b{rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)};
        const Hull hull_a = box_hull(half_a), hull_b = box_hull(half_b);
        const Vec3 pa = rng.in_box(0.5f);
        const Quat qa = rng.rotation(), qb_free = rng.rotation();
        // Half the time b is lined up with a (to within a degree), so faces meet faces.
        const Quat qb = trial % 2 ? qa * Quat::from_axis_angle(rng.direction(), rng.range(0, 0.02f)) : qb_free;
        Convex box_a = Convex::box(half_a, pa, qa), box_b = Convex::box(half_b, rng.in_box(1.2f), qb);
        // Slide b to a shallow overlap.
        const ClosestResult apart = closest(box_a, box_b);
        box_b.pos = box_b.pos - apart.normal * (apart.distance + rng.range(0.002f, 0.03f));
        const Convex as_hull_a = Convex::hull(hull_a, box_a.pos, qa), as_hull_b = Convex::hull(hull_b, box_b.pos, qb);

        Manifold boxes, hulls, mixed;
        const bool hit = collide_convex(box_a, box_b, boxes);
        CHECK(collide_convex(as_hull_a, as_hull_b, hulls) == hit);
        CHECK(collide_convex(box_a, as_hull_b, mixed) == hit);
        if (!hit) continue;
        ++hits;
        patches += boxes.count > 1;
        for (const Manifold* m : {&hulls, &mixed}) {
            CHECK(dot(m->normal, boxes.normal) > 0.9999f);
            CHECK_NEAR(m->depth, boxes.depth, 3e-3);  // (their normals may be a fraction of a degree apart)
            // The same patch. (Up to four of as many as eight candidate points are kept, and which four
            // can differ when the two shapes list their corners in a different order; so compare the
            // patches only when nothing had to be left out.)
            if (boxes.count < 4) CHECK(same_points(boxes, *m, 2e-3f));
            else CHECK(m->count == 4);
        }
    }
    CHECK(hits > 3500);
    CHECK(patches > 1000);
}

// ---- the patch for each kind of hull contact ------------------------------------------------------

TEST(a_tetrahedron_rests_on_three_points) {
    const std::vector<Vec3> corners{{0, 0, 0}, {1, 0, 0}, {0, 0, 1}, {0.3f, 0.8f, 0.3f}};  // one face in the plane y = 0
    const Hull h = built(corners);
    Manifold m;
    CHECK(collide_convex(kSlab, Convex::hull(h, {2, -0.01f, 1}), m));
    CHECK(distance(m.normal, Vec3{0, 1, 0}) < 1e-5f);
    CHECK_NEAR(m.depth, 0.01, 1e-5);
    CHECK(m.count == 3);
    for (const Vec3& c : {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 0, 1}}) {
        const int k = find_point(m, c + Vec3{2, -0.005f, 1}, 1e-4f);
        CHECK(k >= 0);
        if (k >= 0) CHECK_NEAR(m.points[k].depth, 0.01, 1e-5);
    }
}

TEST(a_hexagonal_prism_on_its_end_gets_four_of_its_six_corners) {
    const Hull h = prism(6, 1, 0.5f);
    Manifold m;
    CHECK(collide_convex(kSlab, Convex::hull(h, {0, 0.49f, 0}), m));
    CHECK(distance(m.normal, Vec3{0, 1, 0}) < 1e-5f);
    CHECK_NEAR(m.depth, 0.01, 1e-5);
    CHECK(m.count == 4);
    Vec3 middle;
    for (int i = 0; i < m.count; ++i) {
        // Each is a corner of the hexagon: on the unit circle, at a multiple of 60 degrees.
        const Vec3 p = m.points[i].point;
        CHECK_NEAR(p.y, -0.005, 1e-5);
        CHECK_NEAR(std::sqrt(p.x * p.x + p.z * p.z), 1, 1e-4);
        const Real sixths = std::atan2(p.z, p.x) / (kPi / 3);
        CHECK_NEAR(sixths, std::round(sixths), 1e-3);
        for (int j = 0; j < i; ++j) CHECK(distance(p, m.points[j].point) > 0.5f && m.points[j].id != m.points[i].id);
        middle += p * 0.25f;
    }
    // And the four kept are spread round it: the centre of mass above stays well inside their outline.
    CHECK(std::sqrt(middle.x * middle.x + middle.z * middle.z) < 0.3f);

    // Lying on one of its rectangular sides instead: that rectangle's four corners.
    const Quat on_side = Quat::from_axis_angle({0, 0, 1}, kPi / 2) * Quat::from_axis_angle({0, 1, 0}, kPi / 6);
    const Real to_side = std::cos(kPi / 6);  // from the axis to the middle of a side
    CHECK(collide_convex(kSlab, Convex::hull(h, {0, to_side - 0.01f, 0}, on_side), m));
    CHECK(m.count == 4);
    CHECK_NEAR(m.depth, 0.01, 1e-4);
    for (int i = 0; i < m.count; ++i) {
        CHECK_NEAR(std::fabs(m.points[i].point.x), 0.5, 1e-3);  // the prism's two ends
        CHECK_NEAR(std::fabs(m.points[i].point.z), 0.5, 1e-3);  // the side's two long edges, half a side apart
    }
}

TEST(the_four_points_kept_of_a_many_sided_face_cover_most_of_it) {
    // Four points are all a manifold holds. Of a many-sided face's corners the four kept should be spread
    // right round it, not bunched on one side, or the body would be propped on a sliver of its base. Compare
    // the area they span with the best any four of the corners could do, on faces of several shapes:
    // regular, stretched four to one, and lopsided.
    Lcg rng;
    Real worst = 1;
    for (int shape = 0; shape < 40; ++shape) {
        const int sides = 6 + shape % 10;
        const Real stretch_x = shape % 3 == 1 ? 4.0f : 1.0f, stretch_z = shape % 3 == 2 ? 4.0f : 1.0f;
        std::vector<Vec3> base, points;
        for (int i = 0; i < sides; ++i) {
            // Corners at uneven angles round an ellipse.
            const Real angle = 2 * kPi * (static_cast<Real>(i) + (shape >= 10 ? rng.range(-0.3f, 0.3f) : 0.0f)) / static_cast<Real>(sides);
            base.push_back({stretch_x * std::cos(angle), 0, stretch_z * std::sin(angle)});
            points.push_back(base.back() + Vec3{0, -0.5f, 0});
            points.push_back(base.back() + Vec3{0, 0.5f, 0});
        }
        const Hull h = built(points);
        auto area_of = [](const Vec3* p, int n) {  // of the convex polygon through n points in the xz plane
            Vec3 middle;
            for (int i = 0; i < n; ++i) middle += p[i] / static_cast<Real>(n);
            std::vector<std::pair<Real, int>> round;
            for (int i = 0; i < n; ++i) round.emplace_back(std::atan2(p[i].z - middle.z, p[i].x - middle.x), i);
            std::sort(round.begin(), round.end());
            Real twice = 0;
            for (int i = 0; i < n; ++i) {
                const Vec3 u = p[round[static_cast<size_t>(i)].second] - middle, v = p[round[static_cast<size_t>((i + 1) % n)].second] - middle;
                twice += std::fabs(u.x * v.z - u.z * v.x);
            }
            return twice * 0.5f;
        };
        Real best = 0;
        for (int i = 0; i < sides; ++i)
            for (int j = i + 1; j < sides; ++j)
                for (int k = j + 1; k < sides; ++k)
                    for (int l = k + 1; l < sides; ++l) {
                        const Vec3 four[4] = {base[static_cast<size_t>(i)], base[static_cast<size_t>(j)], base[static_cast<size_t>(k)], base[static_cast<size_t>(l)]};
                        best = std::max(best, area_of(four, 4));
                    }
        for (int trial = 0; trial < 5; ++trial) {
            // Asked both ways round: with the slab first its face is the reference and the choice is made
            // in the slab's axes; with the hull first, in the hull face's own.
            const Convex placed = Convex::hull(h, {rng.range(-1, 1), 0.49f, rng.range(-1, 1)}, Quat::from_axis_angle({0, 1, 0}, rng.range(-3, 3)));
            for (int order = 0; order < 2; ++order) {
                Manifold m;
                CHECK(order ? collide_convex(placed, kSlab, m) : collide_convex(kSlab, placed, m));
                CHECK(m.count == 4);
                if (m.count != 4) continue;
                const Vec3 kept[4] = {m.points[0].point, m.points[1].point, m.points[2].point, m.points[3].point};
                worst = std::min(worst, area_of(kept, 4) / best);
            }
        }
    }
    // The choice is a quick greedy one (deepest, farthest from it, biggest triangle, farthest on the other
    // side), not a search: measured, it never does worse than 65% of the best.
    CHECK(worst > 0.6f);
}

TEST(a_box_on_a_hexagon_is_clipped_to_the_hexagons_outline) {
    // The hexagonal end of a big prism is the reference face this time, and the cube on it overhangs one of
    // the slanted sides: the patch is cut along that side, which no axis-aligned rectangle test would do.
    const Hull h = prism(6, 2, 0.5f);
    const Convex table = Convex::hull(h, {0, -0.5f, 0});  // its top hexagon in the plane y = 0
    // The side between the corners at 0 and 60 degrees is the line x cos30 + z sin30 = 2 cos30.
    const Vec3 out{std::cos(kPi / 6), 0, std::sin(kPi / 6)};
    const Real edge = 2 * std::cos(kPi / 6);
    const Convex cube = Convex::box({0.4f, 0.4f, 0.4f}, out * edge + Vec3{0, 0.39f, 0});  // centred on that side
    Manifold m;
    CHECK(collide_convex(table, cube, m));
    CHECK(distance(m.normal, Vec3{0, 1, 0}) < 1e-5f);
    CHECK(m.count == 4);
    int on_the_cut = 0;
    for (int i = 0; i < m.count; ++i) {
        const Vec3 p = m.points[i].point;
        CHECK(dot(p, out) <= edge + 1e-4f);  // nothing beyond the hexagon's side
        // Nor outside the cube (give or take the 2 mm clip margin, which is taken back square to the
        // hexagon's side and so slightly skew to the cube's).
        CHECK(std::fabs(p.x - out.x * edge) <= 0.4f + 2e-3f && std::fabs(p.z - out.z * edge) <= 0.4f + 2e-3f);
        on_the_cut += std::fabs(dot(p, out) - edge) < 1e-4f;
    }
    CHECK(on_the_cut == 2);  // two of the four are where the cube's outline crosses the side

    // Asked the other way round, the cube's face is the reference and the hexagon is what gets clipped:
    // the same patch.
    Manifold swapped;
    CHECK(collide_convex(cube, table, swapped));
    CHECK(swapped.count == 4);
    for (int i = 0; i < m.count; ++i) CHECK(find_point(swapped, m.points[i].point, 5e-3f) >= 0);
}

TEST(a_capsule_and_a_ball_on_a_hull_face) {
    const Hull h = prism(5, 2, 0.5f);
    const Convex table = Convex::hull(h, {0, -0.5f, 0});
    Manifold m;
    // A capsule lying on the pentagon: a point under each end.
    CHECK(collide_convex(table, Convex::capsule(0.5f, 0.2f, {0.1f, 0.19f, 0.2f}, Quat::from_axis_angle({0, 0, 1}, kPi / 2)), m));
    CHECK(m.count == 2);
    CHECK_NEAR(m.depth, 0.01, 1e-4);
    CHECK(find_point(m, {0.6f, -0.005f, 0.2f}, 1e-3f) >= 0 && find_point(m, {-0.4f, -0.005f, 0.2f}, 1e-3f) >= 0);
    // A ball: one point, under its centre.
    CHECK(collide_convex(table, Convex::sphere(0.3f, {0.5f, 0.28f, -0.4f}), m));
    CHECK(m.count == 1);
    CHECK(distance(m.points[0].point, Vec3{0.5f, -0.01f, -0.4f}) < 1e-3f);
    // A rounded hull (a skin of 0.1 round the prism) is 0.1 taller.
    Convex rounded = table;
    rounded.radius = 0.1f;
    CHECK(collide_convex(rounded, Convex::sphere(0.3f, {0.5f, 0.38f, -0.4f}), m));
    CHECK_NEAR(m.depth, 0.02, 1e-4);
}

TEST(a_cloud_that_was_never_built_into_a_hull_still_gets_one_point) {
    std::vector<Vec3> corners;
    for (int i = 0; i < 8; ++i) corners.push_back({(i & 1) ? 0.5f : -0.5f, (i & 2) ? 0.5f : -0.5f, (i & 4) ? 0.5f : -0.5f});
    Manifold m;
    // On a slab that is itself only a cloud, so that neither shape has a face to offer.
    std::vector<Vec3> slab_corners;
    for (int i = 0; i < 8; ++i) slab_corners.push_back({(i & 1) ? 5.0f : -5.0f, (i & 2) ? 0.0f : -1.0f, (i & 4) ? 5.0f : -5.0f});
    CHECK(collide_convex(Convex::hull(slab_corners.data(), 8, {0, 0, 0}), Convex::hull(corners.data(), 8, {0, 0.49f, 0}), m));
    CHECK(m.count == 1);
    CHECK_NEAR(m.depth, 0.01, 1e-4);
}

// ---- properties of any manifold between hulls ----------------------------------------------------

TEST(every_shallow_manifold_between_hulls_is_inside_both_and_no_deeper_than_the_pair) {
    Lcg rng;
    std::vector<Hull> hulls;
    hulls.push_back(prism(3, 0.6f, 0.4f));
    hulls.push_back(prism(6, 0.5f, 0.3f));
    hulls.push_back(prism(8, 0.6f, 0.2f));
    hulls.push_back(box_hull({0.3f, 0.5f, 0.4f}));
    hulls.push_back(built({{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}}));  // a tetrahedron
    for (int i = 0; i < 5; ++i) hulls.push_back(random_hull(rng));
    int checked = 0, patches = 0, by_count[5] = {};
    for (int trial = 0; trial < 6000; ++trial) {
        const Hull& ha = hulls[static_cast<size_t>(trial) % hulls.size()];
        const Hull& hb = hulls[static_cast<size_t>(trial / 10) % hulls.size()];
        const Convex a = Convex::hull(ha, rng.in_box(0.3f), rng.rotation());
        const int other = trial % 4;  // another hull, a box, a capsule, a ball
        Convex b = other == 0   ? Convex::hull(hb, rng.in_box(1.5f), rng.rotation(), trial % 8 == 0 ? 0.1f : 0.0f)
                   : other == 1 ? Convex::box({rng.range(0.2f, 0.8f), rng.range(0.2f, 0.8f), rng.range(0.2f, 0.8f)}, rng.in_box(1.5f), rng.rotation())
                   : other == 2 ? Convex::capsule(rng.range(0.2f, 0.8f), rng.range(0.1f, 0.3f), rng.in_box(1.5f), rng.rotation())
                                : Convex::sphere(rng.range(0.2f, 0.6f), rng.in_box(1.5f));
        // Half the time, press one of b's features flat against a face of a: turn b so that the contact
        // normal becomes one of a's face normals. (Done by first finding where they would touch.)
        if ((trial / 4) % 2 && other != 3) {
            const Hull::Face& face = ha.faces[static_cast<size_t>(rng.next() * static_cast<Real>(ha.faces.size())) % ha.faces.size()];
            const Vec3 n = rotate(a.q, face.normal);
            const Vec3 from = other == 2 ? rotate(b.q, {1, 0, 0}) : other == 1 ? rotate(b.q, {0, -1, 0}) : rotate(b.q, hb.faces[0].normal);
            // Rotate b so that `from` (one of its own face normals, or a direction across a capsule) points
            // against n: then that face of b lies flat on the face of a.
            const Vec3 axis = cross(from, -n);
            if (axis.length() > 1e-3f) b.q = Quat::from_axis_angle(axis, std::atan2(axis.length(), dot(from, -n))) * b.q;
            b.pos = a.pos + n * 3 + rotate(a.q, ha.point(face, 0)) * 0.3f;
        }
        const ClosestResult apart = closest(a, b);
        b.pos = b.pos - apart.normal * (apart.distance + rng.range(0.001f, 0.02f));

        Manifold m;
        if (!collide_convex(a, b, m)) continue;
        ++checked;
        patches += m.count > 1;
        ++by_count[m.count];
        CHECK(m.count >= 1 && m.count <= Manifold::kMaxPoints);
        CHECK_NEAR(m.normal.length(), 1, 1e-4);
        CHECK_NEAR(m.depth, dot(m.normal, a.support(m.normal) - b.support(-m.normal)), 1e-3);
        for (int i = 0; i < m.count; ++i) {
            const ContactPoint& p = m.points[i];
            CHECK(p.depth >= 0 && p.depth <= m.depth + 1e-4f);
            // A contact point sits halfway through the overlap, below a corner of one shape's face. For a
            // box that is inside the box. For a shape whose sides slope inward from that face (a
            // tetrahedron, most hulls) it can be just outside a sloping side, by less than the overlap.
            CHECK(inside(a, p.point, 3e-3f + m.depth, rng, m.normal));
            CHECK(inside(b, p.point, 3e-3f + m.depth, rng, m.normal));
            for (int j = 0; j < i; ++j) CHECK(m.points[j].id != p.id);
        }
        Convex moved = b;
        moved.pos += m.normal * (m.depth + 2e-3f);
        Manifold none;
        CHECK(!collide_convex(a, moved, none));
    }
    CHECK(checked > 5500);
    CHECK(patches > 1500);                      // plenty were face contacts
    CHECK(by_count[3] > 50 && by_count[4] > 300);  // triangles and bigger faces both
}

TEST(hull_contact_ids_survive_a_small_movement_and_change_with_the_face) {
    Lcg rng;
    const Hull h = prism(6, 0.6f, 0.4f);
    std::vector<std::uint32_t> on_end, on_other_end;
    for (int trial = 0; trial < 100; ++trial) {
        const Quat yaw = Quat::from_axis_angle({0, 1, 0}, rng.range(-3, 3));
        const Vec3 at{rng.range(-2, 2), 0.39f, rng.range(-2, 2)};
        Manifold before, after, flipped;
        CHECK(collide_convex(kSlab, Convex::hull(h, at, yaw), before));
        // Nudged: a millimetre sideways, half a millimetre down, a hair of a turn, and a tilt so slight that
        // it only changes WHICH corner is deepest, by a few hundredths of a millimetre.
        const Quat nudge = Quat::from_axis_angle(rng.direction(), 5e-5f) * Quat::from_axis_angle({0, 1, 0}, 0.002f);
        CHECK(collide_convex(kSlab, Convex::hull(h, at + Vec3{0.001f, -0.0005f, 0.001f}, nudge * yaw), after));
        CHECK(before.count == 4 && after.count == 4);
        for (int i = 0; i < before.count && after.count == before.count; ++i) {
            const int j = find_point(after, before.points[i].point, 5e-3f);
            CHECK(j >= 0);
            if (j >= 0) CHECK(after.points[j].id == before.points[i].id);
        }
        // Turned over onto its other hexagonal end: none of the same ids.
        CHECK(collide_convex(kSlab, Convex::hull(h, at, yaw * Quat::from_axis_angle({1, 0, 0}, kPi)), flipped));
        for (int i = 0; i < before.count; ++i)
            for (int j = 0; j < flipped.count; ++j) CHECK(before.points[i].id != flipped.points[j].id);
    }
}
