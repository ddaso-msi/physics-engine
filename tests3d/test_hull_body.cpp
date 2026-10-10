#include "test.hpp"
#include <phys3d/world.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
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

#define CHECK_MAT(a, b, eps)                                     \
    do {                                                         \
        const Mat3 ma_ = (a), mb_ = (b);                         \
        for (int r_ = 0; r_ < 3; ++r_)                           \
            for (int c_ = 0; c_ < 3; ++c_)                       \
                CHECK_NEAR(ma_.at(r_, c_), mb_.at(r_, c_), eps); \
    } while (0)

std::vector<Vec3> box_corners(Vec3 half, Vec3 shift = {}) {
    std::vector<Vec3> p;
    for (int i = 0; i < 8; ++i) p.push_back(Vec3{(i & 1) ? half.x : -half.x, (i & 2) ? half.y : -half.y, (i & 4) ? half.z : -half.z} + shift);
    return p;
}
std::vector<Vec3> prism_points(int sides, Real r, Real half_height) {
    std::vector<Vec3> p;
    for (int i = 0; i < sides; ++i) {
        const Real a = 2 * kPi * static_cast<Real>(i) / static_cast<Real>(sides);
        p.push_back({r * std::cos(a), -half_height, r * std::sin(a)});
        p.push_back({r * std::cos(a), half_height, r * std::sin(a)});
    }
    return p;
}
std::vector<Vec3> lump(Lcg& rng, Real size = 0.5f, int n = 12) {
    std::vector<Vec3> p;
    for (int i = 0; i < n; ++i) p.push_back(rng.in_box(size));
    return p;
}
// The twelve corners of an icosahedron with edges of length 2 * scale.
std::vector<Vec3> icosahedron(Real scale) {
    const Real g = (1 + std::sqrt(5.0f)) / 2;
    std::vector<Vec3> p;
    for (const Real a : {-1.0f, 1.0f})
        for (const Real b : {-g, g}) {
            p.push_back(Vec3{0, a, b} * scale);
            p.push_back(Vec3{a, b, 0} * scale);
            p.push_back(Vec3{b, 0, a} * scale);
        }
    return p;
}
Hull built(const std::vector<Vec3>& points) {
    Hull h;
    CHECK(build_hull(points.data(), static_cast<int>(points.size()), h));
    return h;
}
Body hull_body(const std::vector<Vec3>& points, Vec3 origin, Quat q = {}, Real density = 1) {
    std::optional<Body> b = Body::solid_hull(points.data(), static_cast<int>(points.size()), density, origin, q);
    CHECK(b.has_value());
    if (!b) return Body::solid_sphere(0.1f, 1, origin);
    b->restitution = 0;
    return *b;
}
bool in_hull(const Hull& h, Vec3 p) {
    for (const Hull::Face& f : h.faces)
        if (dot(f.normal, p) > f.offset) return false;
    return true;
}
// Volume, centre of mass and inertia about that centre of a solid of density 1, the slow way: a grid of
// n^3 equal cells over the hull's bounding box, each inside or not. Independent of the formulas.
MassProperties integrate(const Hull& h, int n) {
    Vec3 lo = h.vertices[0], hi = h.vertices[0];
    for (const Vec3& v : h.vertices)
        for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], v[k]), hi[k] = std::max(hi[k], v[k]);
    const Vec3 cell = (hi - lo) / static_cast<Real>(n);
    std::vector<Vec3> inside;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                const Vec3 p{lo.x + cell.x * (static_cast<Real>(i) + 0.5f), lo.y + cell.y * (static_cast<Real>(j) + 0.5f), lo.z + cell.z * (static_cast<Real>(k) + 0.5f)};
                if (in_hull(h, p)) inside.push_back(p);
            }
    MassProperties out;
    const double each = static_cast<double>(cell.x) * cell.y * cell.z;
    out.volume = static_cast<Real>(each * static_cast<double>(inside.size()));
    double centre[3] = {0, 0, 0};
    for (const Vec3& p : inside)
        for (int k = 0; k < 3; ++k) centre[k] += static_cast<double>(p[k]) / static_cast<double>(inside.size());
    out.centroid = {static_cast<Real>(centre[0]), static_cast<Real>(centre[1]), static_cast<Real>(centre[2])};
    double sum[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (const Vec3& p : inside) {
        const double r[3] = {p.x - centre[0], p.y - centre[1], p.z - centre[2]};
        const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) sum[a][b] += each * ((a == b ? r2 : 0) - r[a] * r[b]);
    }
    auto e = [&](int a, int b) { return static_cast<Real>(sum[a][b]); };
    out.inertia = Mat3{{e(0, 0), e(1, 0), e(2, 0)}, {e(0, 1), e(1, 1), e(2, 1)}, {e(0, 2), e(1, 2), e(2, 2)}};
    return out;
}

World with_floor() {
    World w;
    Body floor = Body::fixed_box({30, 0.5f, 30}, {0, -0.5f, 0});
    floor.restitution = 0;
    w.add(floor);
    return w;
}
void add_room(World& w, Real half) {
    w.add(Body::fixed_box({0.5f, 30, half}, {-half - 0.5f, 30, 0}));
    w.add(Body::fixed_box({0.5f, 30, half}, {half + 0.5f, 30, 0}));
    w.add(Body::fixed_box({half, 30, 0.5f}, {0, 30, -half - 0.5f}));
    w.add(Body::fixed_box({half, 30, 0.5f}, {0, 30, half + 0.5f}));
}
void run(World& w, int steps) {
    for (int i = 0; i < steps; ++i) w.step(kDt);
}

}  // namespace

// ---- weight and inertia of a polyhedron ------------------------------------------------------------

TEST(a_box_shaped_hull_weighs_and_turns_like_the_box) {
    const Vec3 half{0.4f, 0.7f, 1.1f};
    for (const Vec3& shift : {Vec3{0, 0, 0}, Vec3{3, -2, 5}, Vec3{-40, 25, 60}}) {
        const MassProperties m = mass_properties(built(box_corners(half, shift)));
        CHECK_NEAR(m.volume, 8 * 0.4 * 0.7 * 1.1, 1e-4);
        CHECK(distance(m.centroid, shift) < 1e-4f);  // wherever the box is, that is where its centre is...
        CHECK_MAT(m.inertia, box_inertia(m.volume, half), 2e-4);  // ...and its inertia about it is the same
    }
}

TEST(a_tetrahedron_has_a_sixth_of_its_boxs_volume_and_its_centre_at_the_mean_of_its_corners) {
    const std::vector<Vec3> corners{{0.2f, 0.1f, -0.3f}, {1.4f, 0.3f, 0.1f}, {0.5f, 1.6f, 0.2f}, {0.4f, 0.5f, 1.3f}};
    const MassProperties m = mass_properties(built(corners));
    const Vec3 a = corners[1] - corners[0], b = corners[2] - corners[0], c = corners[3] - corners[0];
    CHECK_NEAR(m.volume, std::fabs(dot(a, cross(b, c))) / 6, 1e-5);
    CHECK(distance(m.centroid, (corners[0] + corners[1] + corners[2] + corners[3]) * 0.25f) < 1e-5f);
    // A lopsided solid: turning it about one axis throws it about the others.
    CHECK(std::fabs(m.inertia.at(0, 1)) > 1e-3f);
    CHECK_NEAR(m.inertia.at(0, 1), m.inertia.at(1, 0), 1e-7);
    CHECK_NEAR(m.inertia.at(0, 2), m.inertia.at(2, 0), 1e-7);
    CHECK_NEAR(m.inertia.at(1, 2), m.inertia.at(2, 1), 1e-7);
}

TEST(mass_properties_agree_with_adding_up_little_cells) {
    Lcg rng;
    std::vector<std::vector<Vec3>> shapes{prism_points(5, 0.6f, 0.4f), icosahedron(0.3f), {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    for (int i = 0; i < 6; ++i) shapes.push_back(lump(rng));
    for (std::vector<Vec3> points : shapes) {
        for (Vec3& p : points) p += Vec3{1.5f, -0.5f, 2};  // away from the origin, where a mistake would show
        const Hull h = built(points);
        const MassProperties exact = mass_properties(h), slow = integrate(h, 70);
        CHECK_NEAR(exact.volume, slow.volume, 0.02 * slow.volume);
        CHECK(distance(exact.centroid, slow.centroid) < 0.01f);
        const Real size = slow.inertia.at(0, 0) + slow.inertia.at(1, 1) + slow.inertia.at(2, 2);
        CHECK_MAT(exact.inertia, slow.inertia, 0.02 * size);
    }
}

// ---- the body ------------------------------------------------------------------------------------

TEST(a_hull_body_sits_where_it_was_put_with_its_position_at_its_centre_of_mass) {
    Lcg rng;
    for (int trial = 0; trial < 100; ++trial) {
        std::vector<Vec3> points = lump(rng);
        const Vec3 off_centre = rng.in_box(2);
        for (Vec3& p : points) p += off_centre;  // given in a frame whose origin is nowhere near the solid
        const Vec3 origin = rng.in_box(5);
        const Quat q = rng.rotation();
        const Body b = hull_body(points, origin, q, 3);
        const MassProperties given = mass_properties(built(points));

        CHECK_NEAR(b.mass, 3 * given.volume, 1e-4 * b.mass + 1e-6);
        // Three times as dense: three times as heavy and three times as hard to turn, about the same axes.
        CHECK_MAT(b.inertia_body, given.inertia * 3, 1e-3 * (b.inertia_body.at(0, 0) + b.inertia_body.at(1, 1) + b.inertia_body.at(2, 2)) + 1e-7);
        CHECK(distance(b.pos, origin + rotate(q, given.centroid)) < 1e-4f);
        // The shape is stored about the centre of mass...
        CHECK(mass_properties(*b.shape.hull).centroid.length() < 1e-4f);
        // ...and in the world its corners are exactly where the caller's frame put the given points.
        for (const Vec3& corner : b.shape.hull->vertices) {
            const Vec3 world = apply(b.transform(), corner);
            bool found = false;
            for (const Vec3& p : points) found = found || distance(world, origin + rotate(q, p)) < 1e-4f;
            CHECK(found);
        }
        CHECK_MAT(b.inertia_body * b.inv_inertia_body, Mat3::identity(), 1e-3);
    }
}

TEST(points_with_no_volume_make_no_body_and_a_fixed_hull_cannot_move) {
    const Vec3 flat[4] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    CHECK(!Body::solid_hull(flat, 4, 1, {0, 0, 0}).has_value());
    CHECK(!Body::fixed_hull(flat, 4, {0, 0, 0}).has_value());
    const std::vector<Vec3> corners = box_corners({1, 1, 1});
    const std::optional<Body> wall = Body::fixed_hull(corners.data(), 8, {2, 3, 4});
    CHECK(wall.has_value());
    if (!wall) return;
    CHECK(wall->type == BodyType::Static && wall->inv_mass == 0 && wall->shape.type == Shape::Type::Hull);
    CHECK(wall->inv_inertia_body.at(0, 0) == 0 && wall->inv_inertia_body.at(1, 1) == 0);
    CHECK(distance(wall->pos, Vec3{2, 3, 4}) < 1e-5f);
}

TEST(a_hull_bodys_radii_and_bounding_box_fit_it) {
    Lcg rng;
    for (int trial = 0; trial < 150; ++trial) {
        const Body b = hull_body(trial % 3 == 0 ? prism_points(3 + trial % 6, 0.5f, 0.3f) : lump(rng), rng.in_box(3), rng.rotation());
        const Convex shape = Convex::of(b);
        const AABB box = compute_aabb(b);
        const Vec3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (int k = 0; k < 3; ++k) {
            CHECK_NEAR(box.hi[k], shape.support(axes[k])[k], 1e-4);
            CHECK_NEAR(box.lo[k], shape.support(-axes[k])[k], 1e-4);
        }
        // The inscribed sphere touches the nearest face and pokes through none; the bounding sphere reaches
        // the farthest corner and no further.
        Real nearest_face = 1e9f, farthest_corner = 0;
        for (int k = 0; k < 300; ++k) {
            const Vec3 d = rng.direction();
            nearest_face = std::min(nearest_face, dot(shape.support(d) - b.pos, d));
            farthest_corner = std::max(farthest_corner, distance(shape.support(d), b.pos));
        }
        CHECK(b.inscribed_radius() > 0 && b.inscribed_radius() <= nearest_face + 1e-4f);
        CHECK(b.bounding_radius() >= farthest_corner - 1e-4f);
        Real exact_nearest = 1e9f, exact_farthest = 0;
        for (const Hull::Face& f : b.shape.hull->faces) exact_nearest = std::min(exact_nearest, f.offset);
        for (const Vec3& v : b.shape.hull->vertices) exact_farthest = std::max(exact_farthest, v.length());
        CHECK_NEAR(b.inscribed_radius(), exact_nearest, 1e-6);
        CHECK_NEAR(b.bounding_radius(), exact_farthest, 1e-6);
    }
}

TEST(contains_knows_a_hull_from_its_surroundings) {
    Lcg rng;
    int in = 0, out = 0;
    for (int trial = 0; trial < 3000; ++trial) {
        const Body b = hull_body(lump(rng), rng.in_box(2), rng.rotation());
        const Vec3 p = b.pos + rng.in_box(0.6f);
        const Real gap = closest(Convex::of(b), Convex::sphere(0, p)).distance;  // GJK's opinion
        if (std::fabs(gap) < 2e-3f) continue;
        CHECK(contains(b, p) == (gap < 0));
        (gap < 0 ? in : out) += 1;
    }
    CHECK(in > 300 && out > 1000);
}

// ---- in the world --------------------------------------------------------------------------------

TEST(a_tower_of_box_shaped_hulls_stands_like_a_tower_of_boxes) {
    auto top = [](bool hulls, size_t* awake) {
        World w = with_floor();
        const std::vector<Vec3> corners = box_corners({0.5f, 0.5f, 0.5f});
        for (int i = 0; i < 8; ++i) {
            const Vec3 at{0, 0.5f + static_cast<Real>(i) * 1.001f, 0};
            Body crate = Body::solid_box({0.5f, 0.5f, 0.5f}, 1, at);
            crate.restitution = 0;
            w.add(hulls ? hull_body(corners, at) : crate);
        }
        run(w, 1200);
        *awake = w.stats().awake_bodies;
        return w.bodies[8].pos;
    };
    size_t awake_hulls = 1, awake_boxes = 1;
    const Vec3 as_hulls = top(true, &awake_hulls), as_boxes = top(false, &awake_boxes);
    CHECK_NEAR(as_hulls.y, 7.5, 0.03);
    CHECK(std::fabs(as_hulls.x) < 0.02f && std::fabs(as_hulls.z) < 0.02f);
    CHECK(awake_hulls == 0);
    CHECK(distance(as_hulls, as_boxes) < 0.02f);
}

TEST(a_dropped_tetrahedron_comes_to_rest_on_a_face) {
    // A regular tetrahedron with edges of length a rests with its centre a * sqrt(6) / 12 above the floor.
    const Real a = 1.2f;
    const std::vector<Vec3> corners{Vec3{1, 1, 1} * (a / std::sqrt(8.0f)), Vec3{1, -1, -1} * (a / std::sqrt(8.0f)),
                                    Vec3{-1, 1, -1} * (a / std::sqrt(8.0f)), Vec3{-1, -1, 1} * (a / std::sqrt(8.0f))};
    Lcg rng;
    for (int trial = 0; trial < 6; ++trial) {
        World w = with_floor();
        const int id = w.add(hull_body(corners, {0, 2, 0}, rng.rotation()));
        run(w, 900);
        const Body& b = w.bodies[static_cast<size_t>(id)];
        CHECK_NEAR(b.pos.y, a * std::sqrt(6.0f) / 12, 0.01);
        CHECK(!b.awake);
        CHECK(w.contacts().size() == 1);
        if (w.contacts().size() == 1) CHECK(w.contacts()[0].manifold.count == 3);
    }
}

TEST(a_dropped_twenty_sided_die_comes_to_rest_on_a_face) {
    // An icosahedron with edges of length 2 s has its faces g^2 / sqrt(3) * s from its centre (g the golden
    // ratio) and its corners sqrt(1 + g^2) * s. Resting on a face its centre is at the first height; balanced
    // on an edge or a corner it would be higher.
    const Real s = 0.3f, g = (1 + std::sqrt(5.0f)) / 2;
    Lcg rng;
    for (int trial = 0; trial < 6; ++trial) {
        World w = with_floor();
        Body die = hull_body(icosahedron(s), {0, 2, 0}, rng.rotation());
        die.w = rng.in_box(3);
        die.vel = {rng.range(-1, 1), 0, rng.range(-1, 1)};
        const int id = w.add(die);
        run(w, 1800);
        const Body& b = w.bodies[static_cast<size_t>(id)];
        CHECK_NEAR(b.pos.y, g * g / std::sqrt(3.0f) * s, 0.01);
        CHECK(b.vel.length() < 0.05f);
        if (!w.contacts().empty()) CHECK(w.contacts()[0].manifold.count == 3);
    }
}

TEST(a_prism_stays_standing_on_its_end_and_lying_on_its_side) {
    World w = with_floor();
    const std::vector<Vec3> hexagon = prism_points(6, 0.5f, 0.6f);
    const int standing = w.add(hull_body(hexagon, {-2, 0.61f, 0}));
    // Lying with one rectangular side down: that side is cos(30 degrees) * r from the axis.
    const Quat on_side = Quat::from_axis_angle({0, 0, 1}, kPi / 2) * Quat::from_axis_angle({0, 1, 0}, kPi / 6);
    const int lying = w.add(hull_body(hexagon, {2, 0.5f * std::cos(kPi / 6) + 0.01f, 0}, on_side));
    const Vec3 start = w.bodies[static_cast<size_t>(lying)].pos;
    run(w, 900);
    const Body& up = w.bodies[static_cast<size_t>(standing)];
    const Body& flat = w.bodies[static_cast<size_t>(lying)];
    CHECK_NEAR(up.pos.y, 0.6, 0.01);
    CHECK(rotate(up.q, {0, 1, 0}).y > 0.9999f);  // still upright
    CHECK_NEAR(flat.pos.y, 0.5 * std::cos(kPi / 6), 0.01);
    CHECK(std::fabs(flat.pos.x - start.x) < 0.01f && std::fabs(flat.pos.z - start.z) < 0.01f);  // it has flat sides: it does not roll
    CHECK(w.stats().awake_bodies == 0);
}

TEST(a_pile_of_lumps_and_prisms_settles) {
    World w = with_floor();
    add_room(w, 3);
    const size_t first = w.bodies.size();
    Lcg rng;
    for (int i = 0; i < 36; ++i) {
        const Vec3 at{rng.range(-2, 2), 1.5f + static_cast<Real>(i) * 0.7f, rng.range(-2, 2)};
        w.add(hull_body(i % 3 == 0 ? prism_points(3 + i % 6, 0.45f, 0.3f) : i % 3 == 1 ? lump(rng, 0.45f) : icosahedron(0.25f), at, rng.rotation()));
    }
    Real deepest = 0;
    for (int s = 0; s < 2400; ++s) {
        w.step(kDt);
        if (s > 1200)
            for (const ContactPair& c : w.contacts())
                if (!c.dormant) deepest = std::max(deepest, c.manifold.depth);
    }
    Real fastest = 0;
    for (size_t i = first; i < w.bodies.size(); ++i) {
        const Body& b = w.bodies[i];
        CHECK(b.pos.y > b.inscribed_radius() - 0.05f);  // nothing sank into the floor
        CHECK(std::fabs(b.pos.x) < 3 && std::fabs(b.pos.z) < 3);
        fastest = std::max(fastest, b.vel.length());
    }
    CHECK(fastest < 0.3f);
    CHECK(deepest < 0.05f);
}

TEST(collisions_between_hulls_and_everything_else_conserve_momentum) {
    Lcg rng;
    int collided = 0;
    for (int trial = 0; trial < 120; ++trial) {
        World w;
        w.gravity = {};
        w.allow_sleep = false;
        Body a = hull_body(lump(rng, 0.4f), {-1.5f, 0, 0}, rng.rotation());
        Body b = trial % 4 == 0   ? Body::solid_sphere(rng.range(0.25f, 0.5f), 1, {1.5f, 0, 0})
                 : trial % 4 == 1 ? Body::solid_box({rng.range(0.2f, 0.5f), rng.range(0.2f, 0.5f), rng.range(0.2f, 0.5f)}, 1, {1.5f, 0, 0}, rng.rotation())
                 : trial % 4 == 2 ? Body::solid_capsule(rng.range(0.2f, 0.6f), rng.range(0.15f, 0.3f), 1, {1.5f, 0, 0}, rng.rotation())
                                  : hull_body(prism_points(3 + trial % 5, 0.4f, 0.3f), {1.5f, 0, 0}, rng.rotation());
        a.restitution = b.restitution = 0.5f;
        a.vel = Vec3{3, rng.range(-0.5f, 0.5f), rng.range(-0.5f, 0.5f)};
        b.vel = Vec3{-3, rng.range(-0.5f, 0.5f), rng.range(-0.5f, 0.5f)};
        a.w = rng.in_box(2);
        b.w = rng.in_box(2);
        w.add(a);
        w.add(b);
        auto totals = [&](Vec3& linear, Real& energy) {
            linear = {};
            energy = 0;
            for (const Body& body : w.bodies) {
                linear += body.vel * body.mass;
                energy += body.kinetic_energy();
            }
        };
        Vec3 p0, p1;
        Real e0, e1;
        totals(p0, e0);
        bool touched = false;
        for (int s = 0; s < 180; ++s) {
            w.step(kDt);
            touched = touched || !w.contacts().empty();
        }
        totals(p1, e1);
        collided += touched;
        CHECK(distance(p0, p1) < 1e-4f * (w.bodies[0].mass + w.bodies[1].mass) * 3);
        CHECK(e1 <= e0 * 1.001f);
    }
    CHECK(collided > 90);
}

TEST(a_free_hull_keeps_its_centre_on_a_straight_line_and_gains_no_energy) {
    // The body's position is its centre of mass, so with no force on it that point moves in a straight line
    // however the lopsided solid tumbles around it. (Were `pos` any other point of the body, it would wobble.)
    Lcg rng;
    for (int trial = 0; trial < 20; ++trial) {
        World w;
        w.gravity = {};
        w.allow_sleep = false;
        std::vector<Vec3> points = lump(rng);
        for (Vec3& p : points) p += rng.in_box(1);  // an off-centre frame
        Body b = hull_body(points, {0, 0, 0}, rng.rotation());
        b.vel = {1, 0.5f, -0.25f};
        b.w = rng.in_box(4);
        const Vec3 start = b.pos;
        w.add(b);
        const Real e0 = w.bodies[0].kinetic_energy();
        run(w, 480);
        CHECK(distance(w.bodies[0].pos, start + Vec3{1, 0.5f, -0.25f} * 4) < 1e-3f);
        CHECK(w.bodies[0].kinetic_energy() <= e0 * 1.0001f);
        CHECK(w.bodies[0].kinetic_energy() > e0 * 0.9f);
    }
}

TEST(fast_spinning_hulls_do_not_pass_through_a_thin_wall) {
    int through_with = 0, through_without = 0, stuck = 0;
    for (int trial = 0; trial < 150; ++trial) {
        Lcg rng;
        rng.s = static_cast<std::uint32_t>(trial) * 7919u + 3;
        const std::vector<Vec3> points = trial % 2 ? lump(rng, rng.range(0.1f, 0.3f)) : prism_points(3 + trial % 6, rng.range(0.1f, 0.3f), rng.range(0.03f, 0.3f));
        const Real speed = rng.range(40, 180);
        const Vec3 start{3, rng.range(2, 8), rng.range(-3, 3)}, spin = rng.in_box(40), drift{0, rng.range(-10, 10), rng.range(-10, 10)};
        const Quat q = rng.rotation();
        for (int mode = 0; mode < 2; ++mode) {
            World w;
            w.gravity = {};
            w.allow_sleep = false;
            w.continuous = mode == 1;
            w.add(Body::fixed_box({0.02f, 20, 20}, {8, 5, 0}));
            Body bullet = hull_body(points, start, q);
            bullet.vel = Vec3{speed, 0, 0} + drift;
            bullet.w = spin;
            w.add(bullet);
            run(w, 240);
            (mode == 1 ? through_with : through_without) += w.bodies[1].pos.x > 8;
            if (mode == 1) stuck += w.bodies[1].pos.x < 8 && w.bodies[1].vel.x >= 0.5f;
        }
    }
    CHECK(through_with == 0);
    CHECK(stuck == 0);
    CHECK(through_without > 75);  // without the sweep most do: the test is not vacuous
}
