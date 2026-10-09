#include "test.hpp"
#include <phys/collision.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

using namespace phys;

namespace {

constexpr double kEps = 1e-4;

Body box(Real hx, Real hy, Vec2 pos, Real angle = 0) {
    return Body(Shape::make_polygon(Polygon::box(hx, hy)), pos, angle);
}
Body circle(Real r, Vec2 pos) { return Body(Shape::make_circle(r), pos, 0); }

Polygon regular_polygon(int n, Real radius) {
    std::vector<Vec2> pts;
    for (int k = 0; k < n; ++k) {
        Real a = 2 * kPi * static_cast<Real>(k) / static_cast<Real>(n);
        pts.push_back({radius * std::cos(a), radius * std::sin(a)});
    }
    return *Polygon::from_points(pts);
}

// Deterministic pseudo-random numbers (the engine tests must be reproducible).
struct Lcg {
    std::uint32_t s = 12345;
    Real next() {  // uniform in [0, 1)
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
};

Body random_body(Lcg& rng, Vec2 center) {
    const Vec2 pos = center + Vec2{rng.range(-1.5f, 1.5f), rng.range(-1.5f, 1.5f)};
    const Real angle = rng.range(-kPi, kPi);
    switch (static_cast<int>(rng.next() * 4)) {
        case 0: return Body(Shape::make_polygon(Polygon::box(rng.range(0.4f, 1.2f), rng.range(0.4f, 1.2f))), pos, angle);
        case 1: return Body(Shape::make_polygon(regular_polygon(6, rng.range(0.5f, 1.2f))), pos, angle);
        case 2: {
            const Vec2 tri[3] = {{0, 0}, {rng.range(1.0f, 2.0f), 0}, {rng.range(0.0f, 1.0f), rng.range(0.8f, 1.8f)}};
            return Body(Shape::make_polygon(*Polygon::from_points(tri)), pos, angle);
        }
        default: return Body(Shape::make_circle(rng.range(0.4f, 1.0f)), pos, 0);
    }
}

// Contact points sorted by y so tests do not depend on which one the clipper emits first.
std::vector<ContactPoint> sorted_by_y(const Manifold& m) {
    std::vector<ContactPoint> v(m.points, m.points + m.count);
    std::sort(v.begin(), v.end(), [](const ContactPoint& a, const ContactPoint& b) { return a.point.y < b.point.y; });
    return v;
}

}  // namespace

// ---- circle vs circle ---------------------------------------------------------------------------

TEST(circles_overlap) {
    Manifold m;
    CHECK(collide(circle(1.0f, {0, 0}), circle(0.5f, {1.2f, 0}), m));
    CHECK_NEAR(m.normal.x, 1, kEps);
    CHECK_NEAR(m.normal.y, 0, kEps);
    CHECK_NEAR(m.depth, 0.3, kEps);
    CHECK(m.count == 1);
    CHECK_NEAR(m.points[0].point.x, 0.85, kEps);  // halfway through the 0.3 overlap [0.7, 1.0]
    CHECK_NEAR(m.points[0].depth, 0.3, kEps);
}

TEST(circles_separated) {
    Manifold m;
    CHECK(!collide(circle(1.0f, {0, 0}), circle(0.5f, {1.6f, 0}), m));
}

TEST(circles_concentric_do_not_produce_nan) {
    Manifold m;
    CHECK(collide(circle(1.0f, {2, 3}), circle(0.5f, {2, 3}), m));
    CHECK_NEAR(m.normal.length(), 1, kEps);
    CHECK_NEAR(m.depth, 1.5, kEps);
}

// ---- polygon vs polygon -------------------------------------------------------------------------

TEST(boxes_face_to_face) {
    // Two 2x2 boxes overlapping by 0.2 in x. The contact patch is the whole shared face: two points.
    Manifold m;
    CHECK(collide(box(1, 1, {0, 0}), box(1, 1, {1.8f, 0}), m));
    CHECK_NEAR(m.normal.x, 1, kEps);
    CHECK_NEAR(m.normal.y, 0, kEps);
    CHECK_NEAR(m.depth, 0.2, kEps);
    CHECK(m.count == 2);
    if (m.count == 2) {
        auto pts = sorted_by_y(m);
        CHECK_NEAR(pts[0].point.y, -1, kEps);
        CHECK_NEAR(pts[1].point.y, 1, kEps);
        for (auto& p : pts) {
            CHECK_NEAR(p.point.x, 0.9, kEps);  // halfway between x = 0.8 and x = 1.0
            CHECK_NEAR(p.depth, 0.2, kEps);
        }
    }
}

TEST(boxes_partial_overlap_is_clipped) {
    // B sits up and to the right, so only the upper part of B's left face lies over A's right face.
    Manifold m;
    CHECK(collide(box(1, 1, {0, 0}), box(1, 1, {1.8f, 1.5f}), m));
    CHECK_NEAR(m.normal.x, 1, kEps);
    CHECK_NEAR(m.depth, 0.2, kEps);
    CHECK(m.count == 2);
    if (m.count == 2) {
        auto pts = sorted_by_y(m);
        CHECK_NEAR(pts[0].point.y, 0.5, kEps);  // B's own corner
        CHECK_NEAR(pts[1].point.y, 1.0, kEps);  // clipped where A's face ends
    }
}

TEST(boxes_separated_on_each_axis) {
    Manifold m;
    CHECK(!collide(box(1, 1, {0, 0}), box(1, 1, {2.5f, 0}), m));
    CHECK(!collide(box(1, 1, {0, 0}), box(1, 1, {-2.5f, 0}), m));
    CHECK(!collide(box(1, 1, {0, 0}), box(1, 1, {0, 2.5f}), m));
    CHECK(!collide(box(1, 1, {0, 0}), box(1, 1, {0, -2.5f}), m));
    CHECK(!collide(box(1, 1, {0, 0}), box(1, 1, {2.5f, 2.5f}), m));
}

TEST(diamond_corner_into_box_face) {
    // A square turned 45 degrees touches the box with a single corner: one contact point.
    const Real r2 = std::sqrt(2.0f);
    Manifold m;
    CHECK(collide(box(1, 1, {0, 0}), box(1, 1, {r2 + 0.9f, 0}, kPi / 4), m));
    CHECK_NEAR(m.normal.x, 1, kEps);
    CHECK_NEAR(m.normal.y, 0, kEps);
    CHECK_NEAR(m.depth, 0.1, kEps);
    CHECK(m.count == 1);
    CHECK_NEAR(m.points[0].point.x, 0.95, kEps);
    CHECK_NEAR(m.points[0].point.y, 0, kEps);
}

TEST(swapping_a_and_b_flips_the_normal) {
    Body a = box(1, 1, {0, 0}), b = box(1, 1, {1.8f, 0.3f}, 0.0f);
    Manifold ab, ba;
    CHECK(collide(a, b, ab));
    CHECK(collide(b, a, ba));
    CHECK_NEAR(ab.normal.x, -ba.normal.x, kEps);
    CHECK_NEAR(ab.normal.y, -ba.normal.y, kEps);
    CHECK_NEAR(ab.depth, ba.depth, kEps);
    CHECK(ab.count == ba.count);
}

// ---- polygon vs circle --------------------------------------------------------------------------

TEST(circle_against_box_face) {
    Manifold m;
    CHECK(collide(box(1, 1, {0, 0}), circle(1.0f, {1.8f, 0}), m));
    CHECK_NEAR(m.normal.x, 1, kEps);
    CHECK_NEAR(m.normal.y, 0, kEps);
    CHECK_NEAR(m.depth, 0.2, kEps);
    CHECK(m.count == 1);
    CHECK_NEAR(m.points[0].point.x, 0.9, kEps);  // between the face (x = 1) and the circle's edge (x = 0.8)
}

TEST(circle_against_box_corner) {
    Manifold m;
    CHECK(collide(box(1, 1, {0, 0}), circle(1.0f, {1.5f, 1.5f}), m));
    const double d = std::sqrt(0.5);  // distance from the corner (1,1) to the centre
    CHECK_NEAR(m.normal.x, std::sqrt(0.5), kEps);
    CHECK_NEAR(m.normal.y, std::sqrt(0.5), kEps);
    CHECK_NEAR(m.depth, 1.0 - d, kEps);
    CHECK_NEAR(m.points[0].point.x, 1.5 - std::sqrt(0.5) * (1.0 + d) / 2.0, kEps);
}

TEST(circle_near_corner_but_not_touching) {
    // Close enough to pass the face tests (|x|, |y| each < r from the faces) yet clear of the corner.
    Manifold m;
    CHECK(!collide(box(1, 1, {0, 0}), circle(1.0f, {1.8f, 1.8f}), m));
}

TEST(circle_centre_inside_box) {
    Manifold m;
    CHECK(collide(box(1, 1, {0, 0}), circle(1.0f, {0.2f, 0}), m));
    CHECK_NEAR(m.normal.x, 1, kEps);
    CHECK_NEAR(m.depth, 1.8, kEps);  // 0.8 to the face plus the full radius
}

TEST(circle_first_or_second_gives_opposite_normals) {
    Manifold bc, cb;
    Body b = box(1, 1, {0, 0}), c = circle(1.0f, {1.8f, 0});
    CHECK(collide(b, c, bc));
    CHECK(collide(c, b, cb));
    CHECK_NEAR(bc.normal.x, 1, kEps);
    CHECK_NEAR(cb.normal.x, -1, kEps);
    CHECK_NEAR(bc.depth, cb.depth, kEps);
}

TEST(circle_against_rotated_box) {
    // Box rotated 90 degrees: its 2-wide axis now points up. A circle above it must see the top face.
    Manifold m;
    CHECK(collide(box(2, 1, {0, 0}, kPi / 2), circle(0.5f, {0, 2.3f}), m));
    CHECK_NEAR(m.normal.x, 0, kEps);
    CHECK_NEAR(m.normal.y, 1, kEps);
    CHECK_NEAR(m.depth, 0.2, kEps);
}

// ---- properties over many random pairs -----------------------------------------------------------

// Whatever the shapes, the reported (normal, depth) must be the MINIMUM translation that separates
// them: moving B along the normal by a hair more than `depth` makes them disjoint, and moving it by
// noticeably less leaves them overlapping. This pins down both the normal and the depth exactly.
TEST(random_pairs_report_minimum_translation) {
    Lcg rng;
    int collisions = 0;
    for (int trial = 0; trial < 3000; ++trial) {
        Body a = random_body(rng, {0, 0}), b = random_body(rng, {0, 0});
        Manifold m;
        if (!collide(a, b, m)) continue;
        ++collisions;

        CHECK_NEAR(m.normal.length(), 1, 1e-4);
        CHECK(m.count >= 1 && m.count <= 2);
        CHECK(m.depth >= 0);
        for (int i = 0; i < m.count; ++i) CHECK(m.points[i].depth >= -1e-5f && m.points[i].depth <= m.depth + 1e-3f);

        Body past = b;
        past.pos = b.pos + m.normal * (m.depth + 0.002f);
        Manifold unused;
        CHECK(!collide(a, past, unused));

        if (m.depth > 0.02f) {
            Body short_of = b;
            short_of.pos = b.pos + m.normal * (m.depth - 0.01f);
            CHECK(collide(a, short_of, unused));
        }
    }
    CHECK(collisions > 800);  // the test must actually exercise collisions
}

// If a sample point lies inside both bodies, the narrow phase must report a collision.
TEST(random_pairs_agree_with_point_sampling) {
    Lcg rng;
    rng.s = 99;
    int collisions = 0;
    for (int trial = 0; trial < 600; ++trial) {
        Body a = random_body(rng, {0, 0}), b = random_body(rng, {0, 0});
        bool shared_point = false;
        for (Real y = -4.5f; y <= 4.5f && !shared_point; y += 0.1f)
            for (Real x = -4.5f; x <= 4.5f; x += 0.1f)
                if (a.contains({x, y}) && b.contains({x, y})) {
                    shared_point = true;
                    break;
                }
        Manifold m;
        const bool hit = collide(a, b, m);
        if (hit) ++collisions;
        if (shared_point) CHECK(hit);
    }
    CHECK(collisions > 150);
}

// Rotating and translating the whole scene must not change depth, count or (rotated) normal.
TEST(collision_is_invariant_under_rigid_motion) {
    Lcg rng;
    rng.s = 7;
    int checked = 0;
    for (int trial = 0; trial < 1500; ++trial) {
        Body a = random_body(rng, {0, 0}), b = random_body(rng, {0, 0});
        Manifold m1;
        if (!collide(a, b, m1)) continue;

        const Real phi = rng.range(-kPi, kPi);
        const Rot rot(phi);
        const Vec2 shift{rng.range(-5.0f, 5.0f), rng.range(-5.0f, 5.0f)};
        Body a2 = a, b2 = b;
        a2.pos = rotate(rot, a.pos) + shift;
        b2.pos = rotate(rot, b.pos) + shift;
        a2.set_angle(a.angle + phi);
        b2.set_angle(b.angle + phi);

        Manifold m2;
        if (!collide(a2, b2, m2)) {  // allowed only for razor-thin overlaps lost to float rounding
            CHECK(m1.depth < 1e-3f);
            continue;
        }
        ++checked;
        CHECK_NEAR(m2.depth, m1.depth, 2e-3);
        const Vec2 n1 = rotate(rot, m1.normal);
        CHECK_NEAR(dot(n1, m2.normal), 1, 1e-3);
    }
    CHECK(checked > 400);
}

// ---- feature ids (used to match contact points from one frame to the next) -----------------------

TEST(ids_are_distinct_within_a_manifold_and_stable_while_sliding) {
    // A crate slides along a floor, always overlapping by 1 cm: the same two corners stay in contact,
    // so each point must keep its id as the crate moves.
    Body floor(Shape::make_polygon(Polygon::box(10, 0.5f)), {0, 0}, 0, BodyType::Static);
    std::vector<std::uint32_t> first;
    for (int step = 0; step <= 40; ++step) {
        Body crate_body = box(0.5f, 0.5f, {-2.0f + 0.1f * static_cast<Real>(step), 0.99f});
        Manifold m;
        CHECK(collide(floor, crate_body, m));
        CHECK(m.count == 2);
        if (m.count != 2) continue;
        CHECK(m.points[0].id != m.points[1].id);
        std::vector<std::uint32_t> ids = {m.points[0].id, m.points[1].id};
        std::sort(ids.begin(), ids.end());
        if (step == 0) first = ids;
        CHECK(ids == first);
    }
}

TEST(a_different_pair_of_features_gets_different_ids) {
    Body floor(Shape::make_polygon(Polygon::box(10, 0.5f)), {0, 0}, 0, BodyType::Static);
    Manifold flat, quarter_turn;
    CHECK(collide(floor, box(0.5f, 0.5f, {0, 0.99f}), flat));
    // Turned 90 degrees a square looks the same, but a different edge of it is now on the floor.
    CHECK(collide(floor, box(0.5f, 0.5f, {0, 0.99f}, kPi / 2), quarter_turn));
    CHECK(flat.count == 2 && quarter_turn.count == 2);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) CHECK(flat.points[i].id != quarter_turn.points[j].id);
}

TEST(tipping_a_box_onto_a_corner_keeps_that_corners_id) {
    Body floor(Shape::make_polygon(Polygon::box(10, 0.5f)), {0, 0}, 0, BodyType::Static);
    Manifold flat, tipped;
    CHECK(collide(floor, box(0.5f, 0.5f, {0, 0.99f}), flat));
    // Tilted slightly, only the lower corner still touches the floor: it is the same corner, so the
    // same id as one of the two flat contacts, and it can inherit that point's impulse.
    CHECK(collide(floor, box(0.5f, 0.5f, {0, 1.04f}, 0.12f), tipped));
    CHECK(tipped.count == 1);
    bool shared = false;
    for (int i = 0; i < flat.count; ++i) shared = shared || flat.points[i].id == tipped.points[0].id;
    CHECK(shared);
}

TEST(circle_ids_distinguish_faces_from_corners) {
    Manifold face, corner;
    CHECK(collide(box(1, 1, {0, 0}), circle(1.0f, {1.8f, 0}), face));
    CHECK(collide(box(1, 1, {0, 0}), circle(1.0f, {1.5f, 1.5f}), corner));
    CHECK((face.points[0].id & 0x100) == 0);    // a face index
    CHECK((corner.points[0].id & 0x100) != 0);  // a vertex index
    CHECK(face.points[0].id != corner.points[0].id);
}
