#include "test.hpp"
#include <phys3d/world.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace phys3d;

namespace {

constexpr Real kDt = 1.0f / 120.0f;

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

// Is p inside the shape, or within tol of it? (No direction has it beyond the shape's farthest point.)
bool inside(const Convex& shape, Vec3 p, Real tol, Lcg& rng, Vec3 also) {
    if (dot(p - shape.support(also), also) > tol || dot(p - shape.support(-also), -also) > tol) return false;
    for (int i = 0; i < 40; ++i) {
        const Vec3 d = rng.direction();
        if (dot(p - shape.support(d), d) > tol) return false;
    }
    return true;
}

// The point of the manifold nearest to p, or -1 if none is within tol.
int find_point(const Manifold& m, Vec3 p, Real tol = 2e-3f) {
    for (int i = 0; i < m.count; ++i)
        if (distance(m.points[i].point, p) < tol) return i;
    return -1;
}
bool same_points(const Manifold& x, const Manifold& y, Real tol = 2e-3f) {
    if (x.count != y.count) return false;
    for (int i = 0; i < x.count; ++i) {
        const int j = find_point(y, x.points[i].point, tol);
        if (j < 0 || std::fabs(y.points[j].depth - x.points[i].depth) > tol) return false;
    }
    return true;
}

// kind: 0 sphere, 1 box, 2 capsule, 3 rounded box.
Convex random_shape(int kind, Lcg& rng, Vec3 pos) {
    const Quat q = rng.rotation();
    if (kind == 0) return Convex::sphere(rng.range(0.2f, 0.8f), pos);
    if (kind == 2) return Convex::capsule(rng.range(0.2f, 1), rng.range(0.1f, 0.4f), pos, q);
    Convex box = Convex::box({rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)}, pos, q);
    if (kind == 3) box.radius = 0.1f;
    return box;
}

Body crate(Vec3 p, Vec3 half = {0.5f, 0.5f, 0.5f}, Quat q = {}) {
    Body b = Body::solid_box(half, 1, p, q);
    b.restitution = 0;
    return b;
}
World convex_world_with_floor() {
    World w;
    w.narrowphase = NarrowPhase::Convex;
    Body floor = Body::fixed_box({30, 0.5f, 30}, {0, -0.5f, 0});
    floor.restitution = 0;
    w.add(floor);
    return w;
}
void run(World& w, int steps) {
    for (int i = 0; i < steps; ++i) w.step(kDt);
}

}  // namespace

// ---- the patch for each kind of contact ----------------------------------------------------------

TEST(a_box_resting_on_a_box_gets_its_four_corners) {
    const Convex slab = Convex::box({5, 0.5f, 5}, {0, -0.5f, 0});
    const Convex cube = Convex::box({0.5f, 0.5f, 0.5f}, {1, 0.49f, 2}, Quat::from_axis_angle({0, 1, 0}, 0.3f));
    Manifold m;
    CHECK(collide_convex(slab, cube, m));
    CHECK(distance(m.normal, Vec3{0, 1, 0}) < 1e-5f);
    CHECK_NEAR(m.depth, 0.01, 1e-5);
    CHECK(m.count == 4);
    for (int i = 0; i < 4; ++i) {
        const Vec3 corner = Vec3{1, 0, 2} + rotate(cube.q, {(i & 1) ? 0.5f : -0.5f, 0, (i & 2) ? 0.5f : -0.5f});
        const int k = find_point(m, corner + Vec3{0, -0.005f, 0}, 1e-4f);  // halfway through the 1 cm overlap
        CHECK(k >= 0);
        if (k >= 0) CHECK_NEAR(m.points[k].depth, 0.01, 1e-5);
    }
    // The other way round: the same points, the normal reversed.
    Manifold swapped;
    CHECK(collide_convex(cube, slab, swapped));
    CHECK(distance(swapped.normal, Vec3{0, -1, 0}) < 1e-5f);
    CHECK(same_points(m, swapped, 4e-3f));  // (to within the clip margin: the outline used is the other box's)
}

TEST(a_small_box_overhanging_an_edge_is_clipped_to_the_face_under_it) {
    // A cube on a table, half off the edge at x = 1: the patch is the half that is over the table.
    const Convex table = Convex::box({1, 0.5f, 1}, {0, -0.5f, 0});
    const Convex cube = Convex::box({0.5f, 0.5f, 0.5f}, {1, 0.49f, 0});
    Manifold m;
    CHECK(collide_convex(table, cube, m));
    CHECK(m.count == 4);
    for (const Vec3& corner : {Vec3{0.5f, -0.005f, -0.5f}, Vec3{0.5f, -0.005f, 0.5f}, Vec3{1, -0.005f, -0.5f}, Vec3{1, -0.005f, 0.5f}})
        CHECK(find_point(m, corner, 3e-3f) >= 0);  // (the clip planes sit 2 mm outside the face)
}

TEST(a_capsule_lying_on_a_box_gets_a_point_under_each_end) {
    const Convex slab = Convex::box({5, 0.5f, 5}, {0, -0.5f, 0});
    const Quat along_x = Quat::from_axis_angle({0, 0, 1}, -kPi / 2);  // its own y axis laid along world x
    const Convex capsule = Convex::capsule(1, 0.25f, {0.3f, 0.24f, -0.4f}, along_x);
    Manifold m;
    CHECK(collide_convex(slab, capsule, m));
    CHECK(distance(m.normal, Vec3{0, 1, 0}) < 1e-4f);
    CHECK_NEAR(m.depth, 0.01, 1e-4);
    CHECK(m.count == 2);
    CHECK(find_point(m, {1.3f, -0.005f, -0.4f}, 1e-3f) >= 0);
    CHECK(find_point(m, {-0.7f, -0.005f, -0.4f}, 1e-3f) >= 0);
    if (m.count == 2) CHECK(m.points[0].id != m.points[1].id);
}

TEST(a_tilted_or_upright_capsule_touches_at_its_lower_end_only) {
    const Convex slab = Convex::box({5, 0.5f, 5}, {0, -0.5f, 0});
    // Tilted 10 degrees from lying flat: one end is 0.35 higher than the other.
    const Quat tilted = Quat::from_axis_angle({0, 0, 1}, -kPi / 2 + 0.1745f);
    const Vec3 axis = rotate(tilted, {0, 1, 0});
    const Vec3 low_end = Vec3{0, 0.24f, 0};
    const Convex capsule = Convex::capsule(1, 0.25f, low_end - axis * (axis.y < 0 ? 1.0f : -1.0f), tilted);
    Manifold m;
    CHECK(collide_convex(slab, capsule, m));
    CHECK(m.count == 1);
    CHECK(distance(m.normal, Vec3{0, 1, 0}) < 1e-3f);
    CHECK_NEAR(m.depth, 0.01, 1e-3);
    CHECK(distance(m.points[0].point, Vec3{0, -0.005f, 0}) < 2e-3f);

    const Convex upright = Convex::capsule(1, 0.25f, {2, 1.24f, 1});
    CHECK(collide_convex(slab, upright, m));
    CHECK(m.count == 1);
    CHECK(distance(m.points[0].point, Vec3{2, -0.005f, 1}) < 1e-3f);
}

TEST(a_capsule_hanging_over_an_edge_is_clipped_at_the_edge) {
    const Convex table = Convex::box({1, 0.5f, 1}, {0, -0.5f, 0});
    const Convex capsule = Convex::capsule(1, 0.25f, {0.8f, 0.24f, 0}, Quat::from_axis_angle({0, 0, 1}, -kPi / 2));  // from x = -0.2 to 1.8
    Manifold m;
    CHECK(collide_convex(table, capsule, m));
    CHECK(m.count == 2);
    CHECK(find_point(m, {-0.2f, -0.005f, 0}, 1e-3f) >= 0);
    CHECK(find_point(m, {1, -0.005f, 0}, 3e-3f) >= 0);  // cut off where the table ends
}

TEST(capsules_side_by_side_get_two_points_and_crossed_ones_get_one) {
    // Parallel, both along y; b is shifted up by 0.5 and is shorter. They run alongside from y = -0.2 to 1.
    const Convex a = Convex::capsule(1, 0.3f, {0, 0, 0});
    const Convex b = Convex::capsule(0.7f, 0.2f, {0.45f, 0.5f, 0});
    Manifold m;
    CHECK(collide_convex(a, b, m));
    CHECK(distance(m.normal, Vec3{1, 0, 0}) < 1e-4f);
    CHECK_NEAR(m.depth, 0.05, 1e-4);
    CHECK(m.count == 2);
    // Midway through the overlap: a's surface is at x = 0.3, b's at 0.25.
    CHECK(find_point(m, {0.275f, -0.2f, 0}, 1e-3f) >= 0);
    CHECK(find_point(m, {0.275f, 1, 0}, 1e-3f) >= 0);
    if (m.count == 2) CHECK(m.points[0].id != m.points[1].id);
    // Swapped, the same two points.
    Manifold swapped;
    CHECK(collide_convex(b, a, swapped));
    CHECK(same_points(m, swapped, 1e-3f));

    // Crossed at right angles: one point, between the two axes.
    const Convex across = Convex::capsule(1, 0.2f, {0.45f, 0.3f, 0.1f}, Quat::from_axis_angle({1, 0, 0}, kPi / 2));
    CHECK(collide_convex(a, across, m));
    CHECK(m.count == 1);
    CHECK(distance(m.points[0].point, Vec3{0.275f, 0.3f, 0}) < 1e-3f);

    // End to end, in line: parallel, but they do not run alongside each other at all.
    const Convex above = Convex::capsule(0.5f, 0.3f, {0, 2.05f, 0});
    CHECK(collide_convex(a, above, m));
    CHECK(m.count == 1);
    CHECK(distance(m.normal, Vec3{0, 1, 0}) < 1e-4f);
}

TEST(a_rounded_box_rests_on_its_rounded_face) {
    // A box core of half 0.4 with a 0.1 skin is a unit cube with rounded edges; its flat bottom is the
    // core's face moved out by the radius.
    const Convex slab = Convex::box({5, 0.5f, 5}, {0, -0.5f, 0});
    Convex rounded = Convex::box({0.4f, 0.4f, 0.4f}, {0, 0.48f, 0});
    rounded.radius = 0.1f;
    Manifold m;
    CHECK(collide_convex(slab, rounded, m));
    CHECK(m.count == 4);
    CHECK_NEAR(m.depth, 0.02, 1e-4);
    for (const Vec3& corner : {Vec3{0.4f, -0.01f, 0.4f}, Vec3{-0.4f, -0.01f, 0.4f}, Vec3{0.4f, -0.01f, -0.4f}, Vec3{-0.4f, -0.01f, -0.4f}})
        CHECK(find_point(m, corner, 1e-3f) >= 0);

    // And with the rounded one asked first, so that ITS face is the reference: the same depth and points
    // (the slab's face clipped to the rounded box's flat part, give or take the clip margin).
    Manifold swapped;
    CHECK(collide_convex(rounded, slab, swapped));
    CHECK(swapped.count == 4);
    CHECK_NEAR(swapped.depth, 0.02, 1e-4);
    CHECK(distance(swapped.normal, Vec3{0, -1, 0}) < 1e-5f);
    CHECK(same_points(m, swapped, 4e-3f));

    // A capsule lying on the rounded box's top: both radii count.
    const Convex capsule = Convex::capsule(0.3f, 0.2f, {0, 0.48f + 0.5f + 0.2f - 0.015f, 0}, Quat::from_axis_angle({0, 0, 1}, kPi / 2));
    CHECK(collide_convex(rounded, capsule, m));
    CHECK(m.count == 2);
    CHECK_NEAR(m.depth, 0.015, 1e-4);
    if (m.count == 2) CHECK_NEAR(m.points[0].depth, 0.015, 1e-4);
}

TEST(a_box_gets_new_contact_ids_when_it_rests_on_a_different_face) {
    // Otherwise a box that tips over would start its new contact with the impulses of the old one.
    const Convex slab = Convex::box({5, 0.5f, 5}, {0, -0.5f, 0});
    std::vector<std::uint32_t> seen;
    const Quat turns[6] = {Quat{}, Quat::from_axis_angle({1, 0, 0}, kPi / 2), Quat::from_axis_angle({1, 0, 0}, kPi),
                           Quat::from_axis_angle({1, 0, 0}, -kPi / 2), Quat::from_axis_angle({0, 0, 1}, kPi / 2), Quat::from_axis_angle({0, 0, 1}, -kPi / 2)};
    for (const Quat& q : turns) {
        for (int order = 0; order < 2; ++order) {
            const Convex cube = Convex::box({0.5f, 0.5f, 0.5f}, {0, 0.49f, 0}, q);
            Manifold m;
            CHECK(order ? collide_convex(cube, slab, m) : collide_convex(slab, cube, m));
            CHECK(m.count == 4);
            for (int i = 0; i < m.count; ++i) {
                CHECK(std::find(seen.begin(), seen.end(), m.points[i].id) == seen.end());
                seen.push_back(m.points[i].id);
            }
        }
    }
}

// ---- against the code that knows the shapes -------------------------------------------------------

TEST(the_convex_path_agrees_with_the_box_code) {
    Lcg rng;
    int hits = 0, same_normal = 0, same_patch = 0;
    for (int trial = 0; trial < 8000; ++trial) {
        const Body a = crate(rng.in_box(1), {rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)}, rng.rotation());
        const Body b = crate(rng.in_box(1), {rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)}, rng.rotation());
        Manifold exact, general;
        const bool hit = collide(a, b, exact);
        CHECK(collide_convex(Convex::of(a), Convex::of(b), general) == hit);
        if (!hit) continue;
        ++hits;
        // The box code prefers a face axis unless an edge is clearly shallower; EPA takes the true minimum.
        // Where the two choose the same direction, everything else must match.
        if (dot(exact.normal, general.normal) < 0.9999f) {
            CHECK(general.depth <= exact.depth + 1e-3f);  // the true minimum is never the deeper one
            continue;
        }
        ++same_normal;
        same_patch += same_points(exact, general);
    }
    CHECK(hits > 3000);
    CHECK(same_normal > hits * 7 / 10);
    CHECK(same_patch >= same_normal - same_normal / 500);  // all but a handful (a corner within the 2 mm clip margin)
}

TEST(the_convex_path_agrees_with_the_sphere_code) {
    Lcg rng;
    int hits = 0;
    for (int trial = 0; trial < 3000; ++trial) {
        const bool two_balls = trial % 2 == 0;
        const Body a = two_balls ? Body::solid_sphere(rng.range(0.2f, 0.8f), 1, rng.in_box(0.5f))
                                 : crate(rng.in_box(0.5f), {rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)}, rng.rotation());
        const Body ball = Body::solid_sphere(rng.range(0.2f, 0.8f), 1, rng.in_box(1.5f));
        for (int order = 0; order < 2; ++order) {
            const Body& first = order ? ball : a;
            const Body& second = order ? a : ball;
            Manifold exact, general;
            const bool hit = collide(first, second, exact);
            CHECK(collide_convex(Convex::of(first), Convex::of(second), general) == hit || std::fabs(exact.depth) < 1e-4f);
            if (!hit || dot(exact.normal, general.normal) < 0.999f) continue;  // (a centre deep inside a box can tie)
            ++hits;
            CHECK(general.count == 1);
            CHECK_NEAR(general.depth, exact.depth, 1e-3);
            CHECK(distance(general.points[0].point, exact.points[0].point) < 5e-3f);  // (the normals may differ by a fraction of a degree)
        }
    }
    CHECK(hits > 1200);
}

// ---- properties of any manifold ------------------------------------------------------------------

TEST(every_shallow_manifold_is_inside_both_shapes_and_no_deeper_than_the_pair) {
    // Shallow, because that is what a simulation produces, and because "a contact point lies between the
    // two surfaces" only puts it inside both shapes when the overlap is thinner than the shapes are.
    Lcg rng;
    int checked = 0, patches = 0;
    for (int trial = 0; trial < 6000; ++trial) {
        const Convex a = random_shape(trial % 4, rng, rng.in_box(0.4f));
        Convex b = random_shape((trial / 4) % 4, rng, rng.in_box(1.5f));
        // Half the time, line b up with a (give or take a degree or two), so that faces meet faces, or
        // very nearly do.
        if ((trial / 16) % 2) b.q = a.q * Quat::from_axis_angle(rng.direction(), rng.range(0, 0.03f));
        const ClosestResult apart = closest(a, b);
        b.pos = b.pos - apart.normal * (apart.distance + rng.range(0.001f, 0.03f));  // slid to overlap by 1 to 30 mm
        Manifold m, swapped;
        // (Sliding in along the old normal does not guarantee that much overlap: a shorter way out may
        // open up in another direction. And an overlap under about half a millimetre is within what
        // counts as merely touching, where the two orders of asking may disagree.)
        const bool hit = collide_convex(a, b, m), hit_swapped = collide_convex(b, a, swapped);
        CHECK(hit == hit_swapped || std::max(m.depth, swapped.depth) < 1e-3f);
        if (!hit || !hit_swapped) continue;
        ++checked;
        patches += m.count > 1;
        CHECK(m.count >= 1 && m.count <= Manifold::kMaxPoints);
        CHECK_NEAR(m.normal.length(), 1, 1e-4);
        // The depth goes with the normal: it is exactly how far the two overlap along that direction, by
        // the support functions. (Not EPA's depth along EPA's slightly different normal.)
        CHECK_NEAR(m.depth, dot(m.normal, a.support(m.normal) - b.support(-m.normal)), 1e-3);
        for (int i = 0; i < m.count; ++i) {
            const ContactPoint& p = m.points[i];
            CHECK(p.depth >= 0 && p.depth <= m.depth + 1e-4f);
            // A contact point sits halfway through the overlap there, so it is inside both shapes.
            CHECK(inside(a, p.point, 3e-3f, rng, m.normal));
            CHECK(inside(b, p.point, 3e-3f, rng, m.normal));
            for (int j = 0; j < i; ++j) CHECK(m.points[j].id != p.id);
        }
        // Moving b out by the depth separates the pair (to within the touching tolerance).
        Convex moved = b;
        moved.pos += m.normal * (m.depth + 2e-3f);
        Manifold none;
        CHECK(!collide_convex(a, moved, none));
        // Asked the other way round: the opposite normal, and a depth that goes with it.
        // (Not exactly opposite: when both shapes present a face, each call takes its own first shape's
        // as the reference, and the two faces may be a few degrees out of parallel.)
        CHECK(dot(swapped.normal, m.normal) < -0.99f);
        CHECK_NEAR(swapped.depth, dot(swapped.normal, b.support(swapped.normal) - a.support(-swapped.normal)), 1e-3);
    }
    CHECK(checked > 5000);
    CHECK(patches > 500);  // a fair number were face contacts
}

TEST(contact_ids_survive_a_small_movement) {
    // The solver matches points between frames by id to reuse their impulses, so a resting body that
    // shifts by a millimetre must keep the same ids on the same corners.
    Lcg rng;
    const Convex slab = Convex::box({5, 0.5f, 5}, {0, -0.5f, 0});
    for (int trial = 0; trial < 200; ++trial) {
        const Quat yaw = Quat::from_axis_angle({0, 1, 0}, rng.range(-3, 3));
        const Vec3 at{rng.range(-2, 2), 0.49f, rng.range(-2, 2)};
        const bool capsule = trial % 2 == 1;
        auto shape = [&](Vec3 p, Quat q) {
            return capsule ? Convex::capsule(0.8f, 0.5f, p, q * Quat::from_axis_angle({0, 0, 1}, kPi / 2)) : Convex::box({0.5f, 0.5f, 0.5f}, p, q);
        };
        Manifold before, after;
        CHECK(collide_convex(slab, shape(at, yaw), before));
        const Quat nudge = Quat::from_axis_angle({0, 1, 0}, 0.002f) * yaw;
        CHECK(collide_convex(slab, shape(at + Vec3{0.001f, -0.0005f, 0.001f}, nudge), after));
        CHECK(before.count == (capsule ? 2 : 4));
        CHECK(after.count == before.count);
        for (int i = 0; i < before.count && after.count == before.count; ++i) {
            const int j = find_point(after, before.points[i].point, 5e-3f);
            CHECK(j >= 0);
            if (j >= 0) CHECK(after.points[j].id == before.points[i].id);
        }
    }
}

// ---- in the world --------------------------------------------------------------------------------

TEST(the_world_uses_whichever_narrow_phase_it_is_told_to) {
    // Find a pair of boxes on which the two disagree (the box code prefers a face axis; EPA takes the true
    // shortest way out), and check that each setting of the switch produces its own answer.
    Lcg rng;
    int found = 0;
    for (int trial = 0; trial < 2000 && found < 20; ++trial) {
        const Body a = crate(rng.in_box(1), {rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)}, rng.rotation());
        const Body b = crate(rng.in_box(1), {rng.range(0.2f, 1), rng.range(0.2f, 1), rng.range(0.2f, 1)}, rng.rotation());
        Manifold exact, general;
        if (!collide(a, b, exact) || !collide_convex(Convex::of(a), Convex::of(b), general)) continue;
        if (dot(exact.normal, general.normal) > 0.99f) continue;
        ++found;
        for (NarrowPhase path : {NarrowPhase::Specialised, NarrowPhase::Convex}) {
            World w;
            w.gravity = {};
            w.narrowphase = path;
            w.add(a);
            w.add(b);
            w.step(kDt);
            CHECK(w.contacts().size() == 1);
            if (w.contacts().size() != 1) continue;
            const Vec3 expected = path == NarrowPhase::Convex ? general.normal : exact.normal;
            CHECK(distance(w.contacts()[0].manifold.normal, expected) < 1e-5f);
        }
    }
    CHECK(found == 20);
}

TEST(a_tower_stands_on_the_convex_path_exactly_as_on_the_box_code) {
    auto top_of_tower = [](NarrowPhase path, size_t* awake) {
        World w = convex_world_with_floor();
        w.narrowphase = path;
        for (int i = 0; i < 8; ++i) w.add(crate({0, 0.5f + static_cast<Real>(i) * 1.001f, 0}));
        run(w, 1200);
        *awake = w.stats().awake_bodies;
        return w.bodies[8].pos;
    };
    size_t awake_convex = 1, awake_exact = 1;
    const Vec3 convex = top_of_tower(NarrowPhase::Convex, &awake_convex), exact = top_of_tower(NarrowPhase::Specialised, &awake_exact);
    CHECK_NEAR(convex.y, 7.5, 0.03);
    CHECK(std::fabs(convex.x) < 0.01f && std::fabs(convex.z) < 0.01f);
    CHECK(awake_convex == 0);  // it settles and sleeps
    CHECK(distance(convex, exact) < 1e-3f);
}

TEST(a_pile_of_boxes_and_balls_settles_on_the_convex_path) {
    World w = convex_world_with_floor();
    // Walls, so that nothing rolls away.
    w.add(Body::fixed_box({0.5f, 20, 4}, {-4.5f, 20, 0}));
    w.add(Body::fixed_box({0.5f, 20, 4}, {4.5f, 20, 0}));
    w.add(Body::fixed_box({4, 20, 0.5f}, {0, 20, -4.5f}));
    w.add(Body::fixed_box({4, 20, 0.5f}, {0, 20, 4.5f}));
    const size_t first = w.bodies.size();
    Lcg rng;
    for (int i = 0; i < 40; ++i) {
        const Vec3 p{rng.range(-3, 3), 1.5f + static_cast<Real>(i) * 0.6f, rng.range(-3, 3)};
        if (i % 3 == 0) {
            Body ball = Body::solid_sphere(rng.range(0.3f, 0.5f), 1, p);
            ball.restitution = 0;
            w.add(ball);
        } else {
            w.add(crate(p, {rng.range(0.25f, 0.6f), rng.range(0.25f, 0.6f), rng.range(0.25f, 0.6f)}, rng.rotation()));
        }
    }
    Real deepest = 0;
    for (int s = 0; s < 1800; ++s) {
        w.step(kDt);
        if (s > 600)
            for (const ContactPair& c : w.contacts())
                if (!c.dormant) deepest = std::max(deepest, c.manifold.depth);
    }
    Real fastest = 0;
    for (size_t i = first; i < w.bodies.size(); ++i) {
        const Body& b = w.bodies[i];
        CHECK(b.pos.y > 0.2f);  // nothing sank into the floor
        CHECK(std::fabs(b.pos.x) < 4 && std::fabs(b.pos.z) < 4);
        fastest = std::max(fastest, b.vel.length());
    }
    CHECK(fastest < 0.2f);   // at rest, give or take a ball still rolling
    CHECK(deepest < 0.05f);  // and once settled nothing is squeezed deep into anything
}

TEST(a_box_on_a_slope_is_held_by_friction_on_the_convex_path) {
    World w;
    w.narrowphase = NarrowPhase::Convex;
    const Quat tilt = Quat::from_axis_angle({0, 0, 1}, 0.35f);
    Body ramp = Body::fixed_box({6, 0.2f, 4}, {0, 0, 0}, tilt);
    ramp.friction = 0.6f;
    ramp.restitution = 0;
    w.add(ramp);
    const Vec3 n = rotate(tilt, {0, 1, 0});
    Body grippy = crate(n * 0.7f, {0.5f, 0.5f, 0.5f}, tilt);
    grippy.friction = 0.9f;  // sqrt(0.9 x 0.6) = 0.73 > tan(0.35) = 0.365: it stays
    const int held = w.add(grippy);
    Body icy = crate(n * 0.7f + Vec3{0, 0, 2}, {0.5f, 0.5f, 0.5f}, tilt);
    icy.friction = 0.01f;
    const int slid = w.add(icy);
    const Vec3 start = w.bodies[static_cast<size_t>(held)].pos;
    run(w, 240);
    CHECK(distance(w.bodies[static_cast<size_t>(held)].pos, start) < 0.01f);
    CHECK(w.bodies[static_cast<size_t>(slid)].pos.y < start.y - 0.5f);  // the icy one is well on its way down
    CHECK(w.contacts().size() >= 1 && w.contacts()[0].manifold.count == 4);
}

TEST(continuous_collision_works_on_the_convex_path) {
    World w;
    w.narrowphase = NarrowPhase::Convex;
    w.gravity = {};
    w.add(Body::fixed_box({0.02f, 20, 20}, {8, 5, 0}));
    Body bullet = Body::solid_box({0.12f, 0.2f, 0.15f}, 1, {3, 5, 0}, Quat::from_axis_angle({1, 2, 3}, 0.7f));
    bullet.vel = {90, 4, -3};
    bullet.w = {5, -8, 6};
    w.add(bullet);
    run(w, 120);
    CHECK(w.bodies[1].pos.x < 8);
    CHECK(w.bodies[1].vel.x < 0.5f);
}
