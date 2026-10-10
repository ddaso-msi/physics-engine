#include "test.hpp"
#include <phys3d/world.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>

using namespace phys3d;

namespace {

constexpr Real kG = 9.81f;
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

// Lying down: the capsule's own y axis turned to lie along world x, or z.
const Quat kAlongX = Quat::from_axis_angle({0, 0, 1}, -kPi / 2);
const Quat kAlongZ = Quat::from_axis_angle({1, 0, 0}, kPi / 2);

Body capsule(Real half_length, Real radius, Vec3 p, Quat q = {}) {
    Body b = Body::solid_capsule(half_length, radius, 1, p, q);
    b.restitution = 0;
    return b;
}
World with_floor() {
    World w;
    Body floor = Body::fixed_box({30, 0.5f, 30}, {0, -0.5f, 0});
    floor.restitution = 0;
    w.add(floor);
    return w;
}
void add_room(World& w, Real half) {  // four tall walls around the origin
    w.add(Body::fixed_box({0.5f, 30, half}, {-half - 0.5f, 30, 0}));
    w.add(Body::fixed_box({0.5f, 30, half}, {half + 0.5f, 30, 0}));
    w.add(Body::fixed_box({half, 30, 0.5f}, {0, 30, -half - 0.5f}));
    w.add(Body::fixed_box({half, 30, 0.5f}, {0, 30, half + 0.5f}));
}
void run(World& w, int steps) {
    for (int i = 0; i < steps; ++i) w.step(kDt);
}
Vec3 axis_of(const Body& b) { return rotate(b.q, {0, 1, 0}); }

}  // namespace

// ---- the shape -----------------------------------------------------------------------------------

TEST(a_solid_capsule_has_the_mass_of_its_volume) {
    const Body b = Body::solid_capsule(0.6f, 0.25f, 800, {1, 2, 3});
    const double volume = kPi * 0.25 * 0.25 * (1.2 + 4.0 / 3 * 0.25);  // a cylinder and a sphere
    CHECK_NEAR(b.mass, 800 * volume, 1e-3 * 800 * volume);
    CHECK_NEAR(b.inv_mass * b.mass, 1, 1e-6);
    CHECK(b.shape.type == Shape::Type::Capsule);
    CHECK_NEAR(b.inertia_body.at(1, 1) * b.inv_inertia_body.at(1, 1), 1, 1e-5);
    CHECK_NEAR(b.inscribed_radius(), 0.25, 1e-6);
    CHECK_NEAR(b.bounding_radius(), 0.85, 1e-6);
    const Body fixed = Body::fixed_capsule(0.6f, 0.25f, {0, 0, 0});
    CHECK(fixed.type == BodyType::Static && fixed.inv_mass == 0 && fixed.shape.type == Shape::Type::Capsule);
}

TEST(a_capsules_bounding_box_is_tight) {
    Lcg rng;
    for (int trial = 0; trial < 300; ++trial) {
        const Body b = capsule(rng.range(0.1f, 1.2f), rng.range(0.1f, 0.5f), rng.in_box(3), rng.rotation());
        const AABB box = compute_aabb(b);
        const Convex shape = Convex::of(b);
        // Along each world axis the box reaches exactly as far as the shape does.
        const Vec3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (int k = 0; k < 3; ++k) {
            CHECK_NEAR(box.hi[k], shape.support(axes[k])[k], 1e-4);
            CHECK_NEAR(box.lo[k], shape.support(-axes[k])[k], 1e-4);
        }
    }
}

TEST(contains_knows_a_capsule_from_its_surroundings) {
    Lcg rng;
    int in = 0, out = 0;
    for (int trial = 0; trial < 4000; ++trial) {
        const Body b = capsule(rng.range(0.1f, 1), rng.range(0.1f, 0.5f), rng.in_box(1), rng.rotation());
        const Vec3 p = b.pos + rng.in_box(0.9f);
        // The distance from p to the capsule's axis, worked out in world space.
        const Vec3 axis = axis_of(b);
        const Real along = std::clamp(dot(p - b.pos, axis), -b.shape.half_length, b.shape.half_length);
        const Real from_axis = distance(p, b.pos + axis * along);
        if (std::fabs(from_axis - b.shape.radius) < 1e-4f) continue;  // too close to call
        CHECK(contains(b, p) == (from_axis < b.shape.radius));
        (from_axis < b.shape.radius ? in : out) += 1;
    }
    CHECK(in > 200 && out > 1000);
}

TEST(collide_sends_capsules_through_the_convex_path) {
    Lcg rng;
    int hits = 0;
    for (int trial = 0; trial < 1500; ++trial) {
        const Body a = capsule(rng.range(0.2f, 1), rng.range(0.1f, 0.4f), rng.in_box(0.5f), rng.rotation());
        Body other = trial % 3 == 0   ? Body::solid_sphere(rng.range(0.2f, 0.6f), 1, rng.in_box(1.2f))
                     : trial % 3 == 1 ? Body::solid_box({rng.range(0.2f, 0.8f), rng.range(0.2f, 0.8f), rng.range(0.2f, 0.8f)}, 1, rng.in_box(1.2f), rng.rotation())
                                      : capsule(rng.range(0.2f, 1), rng.range(0.1f, 0.4f), rng.in_box(1.2f), rng.rotation());
        for (int order = 0; order < 2; ++order) {
            const Body& first = order ? other : a;
            const Body& second = order ? a : other;
            Manifold got, expected;
            const bool hit = collide(first, second, got);
            CHECK(hit == collide_convex(Convex::of(first), Convex::of(second), expected));
            if (!hit) continue;
            ++hits;
            CHECK(got.count == expected.count && distance(got.normal, expected.normal) == 0 && got.depth == expected.depth);
        }
    }
    CHECK(hits > 500);
}

// ---- in the world --------------------------------------------------------------------------------

TEST(a_capsule_dropped_flat_rests_on_two_points_and_sleeps) {
    World w = with_floor();
    const int log = w.add(capsule(0.8f, 0.3f, {0, 1.5f, 0}, kAlongX));
    run(w, 600);
    const Body& b = w.bodies[static_cast<size_t>(log)];
    CHECK_NEAR(b.pos.y, 0.3, 0.01);
    CHECK(std::fabs(axis_of(b).y) < 1e-3f);  // still lying flat
    CHECK(!b.awake);
    CHECK(w.contacts().size() == 1);
    if (w.contacts().size() == 1) CHECK(w.contacts()[0].manifold.count == 2);
}

TEST(a_capsule_dropped_at_an_angle_ends_up_lying_down) {
    World w = with_floor();
    add_room(w, 4);
    const int first = static_cast<int>(w.bodies.size());
    Lcg rng;
    for (int i = 0; i < 6; ++i)
        w.add(capsule(0.7f, 0.25f, {static_cast<Real>(i % 3) * 2.2f - 2.2f, 1.5f, static_cast<Real>(i / 3) * 2.5f - 1.2f},
                      Quat::from_axis_angle(rng.direction(), rng.range(0.3f, 1.2f))));
    run(w, 1200);
    for (size_t i = static_cast<size_t>(first); i < w.bodies.size(); ++i) {
        const Body& b = w.bodies[i];
        CHECK_NEAR(b.pos.y, 0.25, 0.01);
        CHECK(std::fabs(axis_of(b).y) < 0.01f);
        CHECK(b.vel.length() < 0.5f);  // lying, though it may still be rolling gently
    }
}

TEST(a_capsule_rolls_down_a_slope_at_the_textbook_acceleration) {
    // Rolling without slipping, a = g sin(theta) / (1 + I / (m r^2)), with I the inertia about the axis it
    // rolls on. A capsule is a cylinder and a sphere's worth of ends, so I / (m r^2) lies between 2/5 and
    // 1/2, by how much of its mass is in each.
    const Real theta = 0.2f, half_length = 0.6f, radius = 0.3f;
    World w;
    const Quat tilt = Quat::from_axis_angle({0, 0, 1}, theta);  // the slope falls toward -x
    Body ramp = Body::fixed_box({30, 0.5f, 5}, {0, 0, 0}, tilt);
    ramp.friction = 1;
    ramp.restitution = 0;
    w.add(ramp);
    const Vec3 normal = rotate(tilt, {0, 1, 0}), downhill = -rotate(tilt, {1, 0, 0});
    Body log = capsule(half_length, radius, rotate(tilt, {20, 0, 0}) + normal * (0.5f + radius), tilt * kAlongZ);
    log.friction = 1;
    const int id = w.add(log);
    const Body& b = w.bodies[static_cast<size_t>(id)];
    const Real about_axis = b.inertia_body.at(1, 1) / (b.mass * radius * radius);
    CHECK(about_axis > 0.4f && about_axis < 0.5f);
    const Real expected = kG * std::sin(theta) / (1 + about_axis);

    run(w, 60);
    const Real v0 = dot(b.vel, downhill);
    run(w, 240);
    const Real v1 = dot(b.vel, downhill);
    CHECK_NEAR((v1 - v0) / 2.0, expected, 0.02 * expected);
    // Rolling, not sliding: the surface speed from the spin matches the speed of travel.
    CHECK_NEAR(std::fabs(dot(b.w, axis_of(b))) * radius, v1, 0.01 * v1);
    CHECK(std::fabs(b.pos.z) < 0.05f);  // and it runs straight
}

TEST(collisions_between_capsules_and_everything_else_conserve_momentum) {
    Lcg rng;
    int collided = 0;
    for (int trial = 0; trial < 150; ++trial) {
        World w;
        w.gravity = {};
        w.allow_sleep = false;
        Body a = capsule(rng.range(0.3f, 0.8f), rng.range(0.15f, 0.35f), {-1.5f, 0, 0}, rng.rotation());
        Body b = trial % 3 == 0   ? Body::solid_sphere(rng.range(0.25f, 0.5f), 1, {1.5f, 0, 0})
                 : trial % 3 == 1 ? Body::solid_box({rng.range(0.2f, 0.5f), rng.range(0.2f, 0.5f), rng.range(0.2f, 0.5f)}, 1, {1.5f, 0, 0}, rng.rotation())
                                  : capsule(rng.range(0.3f, 0.8f), rng.range(0.15f, 0.35f), {1.5f, 0, 0}, rng.rotation());
        a.restitution = b.restitution = 0.5f;
        a.vel = Vec3{3, rng.range(-0.5f, 0.5f), rng.range(-0.5f, 0.5f)};
        b.vel = Vec3{-3, rng.range(-0.5f, 0.5f), rng.range(-0.5f, 0.5f)};
        a.w = rng.in_box(2);
        b.w = rng.in_box(2);
        w.add(a);
        w.add(b);
        auto totals = [&](Vec3& linear, Vec3& angular, Real& energy) {
            linear = angular = {};
            energy = 0;
            for (const Body& body : w.bodies) {
                linear += body.vel * body.mass;
                angular += cross(body.pos, body.vel) * body.mass + body.angular_momentum();
                energy += body.kinetic_energy();
            }
        };
        Vec3 p0, l0, p1, l1;
        Real e0, e1;
        totals(p0, l0, e0);
        bool touched = false;
        for (int s = 0; s < 180; ++s) {
            w.step(kDt);
            touched = touched || !w.contacts().empty();
        }
        totals(p1, l1, e1);
        collided += touched;
        const Real scale = (w.bodies[0].mass + w.bodies[1].mass) * 3;
        CHECK(distance(p0, p1) < 1e-4f * scale);
        // Angular momentum only roughly. The contact impulses are equal and opposite at one point, so they
        // cannot change it; what does is the gyroscopic step, once an off-centre hit has set a long thin
        // body tumbling at up to 19 rad/s. Measured: up to 5% here, most of it in free flight afterwards.
        CHECK(distance(l0, l1) < 0.08f * scale);
        CHECK(e1 <= e0 * 1.001f);
    }
    CHECK(collided > 120);
}

TEST(a_pile_of_capsules_settles) {
    World w = with_floor();
    add_room(w, 3);
    const size_t first = w.bodies.size();
    Lcg rng;
    for (int i = 0; i < 30; ++i)
        w.add(capsule(rng.range(0.3f, 0.7f), rng.range(0.2f, 0.3f), {rng.range(-2, 2), 1.5f + static_cast<Real>(i) * 0.7f, rng.range(-2, 2)}, rng.rotation()));
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
        CHECK(b.pos.y > 0.15f);  // nothing sank into the floor
        CHECK(std::fabs(b.pos.x) < 3 && std::fabs(b.pos.z) < 3);
        fastest = std::max(fastest, b.vel.length());
    }
    CHECK(fastest < 0.3f);
    CHECK(deepest < 0.05f);
}

TEST(logs_stacked_crosswise_stay_stacked) {
    // Two logs on the floor, two across them, two across those: every contact between logs is a single
    // point where two round surfaces cross, held by friction alone.
    World w = with_floor();
    const Real r = 0.25f, h = 1.0f;
    for (int layer = 0; layer < 3; ++layer)
        for (const Real side : {-0.6f, 0.6f}) {
            const bool along_x = layer % 2 == 0;
            Body log = capsule(h, r, along_x ? Vec3{0, r + 2 * r * static_cast<Real>(layer), side} : Vec3{side, r + 2 * r * static_cast<Real>(layer), 0},
                               along_x ? kAlongX : kAlongZ);
            log.friction = 0.8f;
            w.add(log);
        }
    run(w, 1200);
    for (int layer = 0; layer < 3; ++layer)
        for (int k = 0; k < 2; ++k) {
            const Body& b = w.bodies[static_cast<size_t>(1 + layer * 2 + k)];
            CHECK_NEAR(b.pos.y, r + 2 * r * static_cast<Real>(layer), 0.02);
            CHECK(std::fabs(axis_of(b).y) < 0.01f);
            CHECK(std::fabs(std::fabs(layer % 2 == 0 ? b.pos.z : b.pos.x) - 0.6f) < 0.02f);  // has not rolled off
        }
    CHECK(w.stats().awake_bodies == 0);
}

TEST(a_fast_spinning_capsule_does_not_pass_through_a_thin_wall) {
    int through_with = 0, through_without = 0;
    for (int trial = 0; trial < 120; ++trial) {
        Lcg rng;
        rng.s = static_cast<std::uint32_t>(trial) * 7919u + 3;
        const Real speed = rng.range(40, 150), radius = rng.range(0.06f, 0.15f), half_length = rng.range(0.1f, 0.5f);
        const Vec3 start{3, rng.range(2, 8), rng.range(-3, 3)}, spin = rng.in_box(20);
        const Quat q = rng.rotation();
        for (int mode = 0; mode < 2; ++mode) {
            World w;
            w.gravity = {};
            w.continuous = mode == 1;
            w.add(Body::fixed_box({0.02f, 20, 20}, {8, 5, 0}));
            Body bullet = capsule(half_length, radius, start, q);
            bullet.vel = {speed, rng.range(-10, 10), rng.range(-10, 10)};
            bullet.w = spin;
            w.add(bullet);
            run(w, 60);
            (mode == 1 ? through_with : through_without) += w.bodies[1].pos.x > 8;
        }
    }
    CHECK(through_with == 0);
    CHECK(through_without > 60);  // and without the sweep most do: the test is not vacuous
}

namespace {

// One of a family of nasty shots: needle-thin to stubby capsules, walls from 1 to 20 cm, up to 200 m/s and
// 80 rad/s about every axis.
World harsh_capsule_shot(int trial) {
    Lcg rng;
    rng.s = static_cast<std::uint32_t>(trial) * 104729u + 11;
    const Real speed = rng.range(20, 200), radius = rng.range(0.02f, 0.15f), half_length = rng.range(0.1f, 0.8f), thick = rng.range(0.01f, 0.2f);
    World w;
    w.gravity = {};
    w.allow_sleep = false;
    w.add(Body::fixed_box({thick * 0.5f, 60, 60}, {8, 5, 0}));
    const Vec3 start{3, rng.range(2, 8), rng.range(-3, 3)};
    const Quat q = Quat::from_axis_angle({rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)}, rng.range(0, 3));
    Body b = Body::solid_capsule(half_length, radius, 1, start, q);
    b.vel = {speed, rng.range(-40, 40), rng.range(-40, 40)};
    b.w = {rng.range(-80, 80), rng.range(-80, 80), rng.range(-80, 80)};
    w.add(b);
    return w;
}
// 0 fine, 1 got through, 2 still this side and still heading in after two seconds, 3 gained energy.
int outcome(World w) {
    const Real e0 = w.bodies[1].kinetic_energy();
    Real peak = 0;
    for (int s = 0; s < 240; ++s) {
        w.step(kDt);
        peak = std::max(peak, w.bodies[1].kinetic_energy());
    }
    const Body& b = w.bodies[1];
    return b.pos.x > 8.3f ? 1 : (b.pos.x < 8 && b.vel.x >= 0.5f) ? 2 : peak > e0 * 1.001f ? 3 : 0;
}

}  // namespace

TEST(harsh_capsule_shots_neither_leak_nor_stick_nor_gain_energy) {
    int through = 0, stuck = 0, gained = 0;
    for (int trial = 0; trial < 500; ++trial) {
        const int result = outcome(harsh_capsule_shot(trial));
        through += result == 1;
        stuck += result == 2;
        gained += result == 3;
    }
    CHECK(through == 0);
    CHECK(stuck == 0);
    CHECK(gained == 0);
}

TEST(capsule_shots_that_once_got_through_now_bounce_off) {
    // Found by stress runs of up to 200000 shots, against earlier versions of the sweep's rule for a body
    // already touching the wall.
    //   1671 to 4147 and 26232 to 50279: a thin capsule arrives end-on and turns side-on as it goes in, so
    //   its farthest point gets no deeper while all of it passes through. The rule watched the farthest
    //   point; it now watches the centre.
    //   21587: the touching end is swinging away from the wall at 340 rad/s, the solver sees a separating
    //   contact and does nothing, and the capsule crept through a quarter radius per step. Being held at the
    //   cap is now treated as the impact it is.
    for (int trial : {1671, 1741, 1894, 3524, 3808, 4147, 26232, 29264, 50279, 21587}) {
        World w = harsh_capsule_shot(trial);
        for (int s = 0; s < 240; ++s) w.step(kDt);
        CHECK(w.bodies[1].pos.x < 8);
        CHECK(w.bodies[1].vel.x < 0.5f);
    }
}

TEST(a_thin_body_spinning_against_a_thin_wall_never_gets_its_centre_into_it) {
    // The body starts already touching the wall, lying against it or leaning on one end, with up to
    // 80 rad/s of spin and little speed toward the wall: it is the spin that would take it through. These
    // particular ones got their centres into or past the wall when the sweep let a body in contact move a
    // full inscribed radius further in per step instead of a quarter.
    for (int trial : {265, 586, 1001, 1321, 2057, 3161, 3479, 4297, 5678, 7046, 8929, 9110, 9290, 10985}) {
        Lcg rng;
        rng.s = static_cast<std::uint32_t>(trial) * 2654435761u + 7;
        const Real radius = rng.range(0.02f, 0.15f), half = rng.range(0.1f, 0.8f), thick = rng.range(0.01f, 0.1f);
        const bool box = trial % 2;
        World w;
        w.gravity = {};
        w.allow_sleep = false;
        w.add(Body::fixed_box({thick * 0.5f, 60, 60}, {8, 5, 0}));
        const Real lean = trial % 4 < 2 ? 0 : rng.range(0.1f, 1.4f);
        const Quat q = Quat::from_axis_angle({0, 0, 1}, lean) * Quat::from_axis_angle({1, 0, 0}, rng.range(0, 3));
        const Vec3 axis = rotate(q, {0, 1, 0});
        const Real reach = box ? std::fabs(axis.x) * half + radius * 1.5f : std::fabs(axis.x) * half + radius;
        const Vec3 start{8 - thick * 0.5f - reach + 0.002f, 5, 0};
        Body b = box ? Body::solid_box({radius, half, radius}, 1, start, q) : Body::solid_capsule(half, radius, 1, start, q);
        b.vel = {rng.range(-2, 10), rng.range(-5, 5), rng.range(-5, 5)};
        b.w = {rng.range(-80, 80), rng.range(-80, 80), rng.range(-80, 80)};
        w.add(b);
        const Real e0 = w.bodies[1].kinetic_energy();
        Real deepest = 0, peak = 0;
        for (int s = 0; s < 240; ++s) {
            w.step(kDt);
            deepest = std::max(deepest, w.bodies[1].pos.x);
            peak = std::max(peak, w.bodies[1].kinetic_energy());
        }
        CHECK(deepest < 8 - thick * 0.5f);
        CHECK(peak <= e0 * 1.001f);
    }
}

TEST(a_body_held_by_the_sweep_bounces_as_the_bouncier_of_the_two_says) {
    // Shot 1802 of the box stress is one the sweep has to hold at the wall and turn round itself. With a
    // dead body and a dead wall it stops; with a dead body and a lively wall it comes back.
    auto speed_back = [](Real wall_restitution) {
        Lcg rng;
        rng.s = 1802u * 104729u + 11;
        const Real speed = rng.range(30, 200), thick = rng.range(0.01f, 0.1f), size = rng.range(0.08f, 0.4f);
        const Vec3 vel{speed, rng.range(-40, 40), rng.range(-40, 40)};
        const Vec3 spin{rng.range(-30, 30), rng.range(-30, 30), rng.range(-30, 30)};
        const Vec3 start{3, rng.range(2, 8), rng.range(-3, 3)};
        const Quat q = Quat::from_axis_angle({rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)}, rng.range(0, 3));
        const Vec3 half{size, size * rng.range(0.3f, 1.6f), size * rng.range(0.3f, 1.6f)};
        World w;
        w.gravity = {};
        Body wall = Body::fixed_box({thick * 0.5f, 60, 60}, {8, 5, 0});
        wall.restitution = wall_restitution;
        w.add(wall);
        Body b = Body::solid_box(half, 1, start, q);
        b.restitution = 0;
        b.vel = vel;
        b.w = spin;
        w.add(b);
        for (int s = 0; s < 240; ++s) w.step(kDt);
        CHECK(w.bodies[1].pos.x < 8);
        return -w.bodies[1].vel.x;
    };
    const Real dead = speed_back(0), lively = speed_back(0.8f);
    CHECK(std::fabs(dead) < 1e-3f);  // stopped, to rounding
    CHECK(lively > dead + 1);
}

TEST(a_capsule_on_a_hinge_swings_like_any_other_body) {
    // A capsule hung by one end is a compound pendulum: T = 2 pi sqrt(I_pivot / (m g d)), with d from the
    // pivot to the centre and I_pivot the transverse inertia moved out by d.
    const Real half_length = 0.6f, radius = 0.1f, start = 0.1f;
    World w;
    const Real d = half_length + radius;  // pivot at the very tip
    const int rod = w.add(capsule(half_length, radius, {d * std::sin(start), -d * std::cos(start), 0}, Quat::from_axis_angle({0, 0, 1}, start)));
    w.add_joint(Joint::hinge(w.bodies, -1, rod, {0, 0, 0}, {0, 0, 1}));
    const Body& b = w.bodies[static_cast<size_t>(rod)];
    const double about_pivot = b.inertia_body.at(2, 2) + b.mass * d * d;
    const double expected = 2 * kPi * std::sqrt(about_pivot / (b.mass * kG * d));
    double first = -1, last = -1;
    int swings = 0;
    Real before = b.pos.x;
    const Real dt = 1.0f / 480.0f;
    for (int i = 1; i <= 480 * 8; ++i) {
        w.step(dt);
        if (before < 0 && b.pos.x >= 0) {
            const double t = (i - 1 + before / (before - b.pos.x)) * static_cast<double>(dt);
            if (first < 0) first = t;
            else ++swings;
            last = t;
        }
        before = b.pos.x;
    }
    CHECK(swings >= 3);
    if (swings >= 3) CHECK_NEAR((last - first) / swings, expected, 0.005 * expected);
}
