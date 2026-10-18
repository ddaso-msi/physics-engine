// Rays, objectives, and the level life cycle, exercised through real levels.
#include "test.hpp"
#include <levels/knock_down.hpp>
#include <levels/levels.hpp>
#include <levels/pendulum.hpp>
#include <levels/tower.hpp>
#include <objectives.hpp>
#include <raycast.hpp>

#include <cmath>
#include <cstring>

using namespace phys3d;
using namespace puzzle3d;

namespace {
void run(Level& l, int steps) {
    for (int i = 0; i < steps; ++i) l.step();
}
Ray ray(Vec3 from, Vec3 to) { return {from, (to - from).normalized()}; }
}  // namespace

// ---- rays ----

TEST(ray_hits_a_sphere_where_it_enters) {
    const Body b = Body::solid_sphere(1, 1, {0, 0, -5});
    const auto h = raycast(b, ray({0, 0, 0}, {0, 0, -5}));
    CHECK(h.has_value());
    CHECK_NEAR(h->t, 4.0, 1e-4);
    CHECK_NEAR(h->normal.z, 1.0, 1e-4);
    CHECK(!raycast(b, ray({0, 2, 0}, {0, 2, -5})).has_value());  // passes over it
    CHECK(!raycast(b, ray({0, 0, 0}, {0, 0, 5})).has_value());   // behind the start
    const auto inside_start = raycast(b, ray({0, 0, -5}, {0, 0, -9}));
    CHECK(inside_start && inside_start->t == 0);
}

TEST(ray_hits_a_turned_box_on_its_face) {
    const Quat q = Quat::from_axis_angle({0, 1, 0}, 0.5f);
    const Body b = Body::solid_box({1, 0.5f, 2}, 1, {3, 1, 0}, q);
    // Straight down onto its top: the top is at y = 1.5 whatever the turn about y.
    const auto h = raycast(b, ray({3, 10, 0}, {3, 0, 0}));
    CHECK(h.has_value());
    CHECK_NEAR(h->point.y, 1.5, 1e-4);
    CHECK_NEAR(h->normal.y, 1.0, 1e-4);
    // Along the box's own long axis, from outside: enters at its end, 2 m from the centre.
    const Vec3 axis = rotate(q, Vec3{0, 0, 1});
    const auto e = raycast(b, {b.pos + axis * 6, -axis});
    CHECK(e.has_value());
    CHECK_NEAR(e->t, 4.0, 1e-3);
    // Just past a corner it misses.
    const Vec3 side = rotate(q, Vec3{1, 0, 0});
    CHECK(!raycast(b, {b.pos + axis * 6 + side * 1.05f, -axis}).has_value());
}

TEST(ray_hits_a_capsule_on_its_side_and_its_end) {
    const Body c = Body::solid_capsule(1, 0.5f, 1, {0, 0, 0});  // axis along y, ends at y = +-1.5
    const auto side = raycast(c, ray({5, 0.5f, 0}, {0, 0.5f, 0}));
    CHECK(side && std::fabs(side->t - 4.5f) < 1e-3f);
    const auto top = raycast(c, ray({0, 5, 0}, {0, 0, 0}));
    CHECK(top && std::fabs(top->point.y - 1.5f) < 1e-3f);
    CHECK(!raycast(c, ray({5, 3, 0}, {0, 3, 0})).has_value());
}

TEST(ray_hits_a_hull) {
    const Vec3 pts[] = {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {0, 0, 2}};  // a corner of a cube, cut off
    const std::optional<Body> h = Body::solid_hull(pts, 4, 1, {0, 0, 0});
    CHECK(h.has_value());
    const auto hit = raycast(*h, ray({0.3f, 0.3f, 10}, {0.3f, 0.3f, 0}));
    CHECK(hit && std::fabs(hit->point.z - 1.4f) < 1e-3f);  // the slanted face x + y + z = 2
    CHECK(!raycast(*h, ray({1.5f, 1.5f, 10}, {1.5f, 1.5f, 0})).has_value());
}

TEST(world_ray_finds_the_nearest_and_can_skip_one) {
    World w;
    add_ground(w);
    const int near = w.add(Body::solid_box({0.5f, 0.5f, 0.5f}, 1, {0, 0.5f, 0}));
    const int far = w.add(Body::solid_box({0.5f, 0.5f, 0.5f}, 1, {0, 0.5f, -3}));
    const Ray r = ray({0, 0.5f, 10}, {0, 0.5f, 0});
    CHECK(raycast(w, r)->body == near);
    CHECK(raycast(w, r, near)->body == far);
    const auto down = raycast(w, ray({5, 5, 5}, {5, 0, 5}));
    CHECK(down && down->body == 0 && std::fabs(down->point.y) < 1e-4f);  // the ground
    CHECK_NEAR(hit_plane_y(ray({0, 5, 0}, {1, 4, 0}), 0)->x, 5.0, 1e-4);
    CHECK(!hit_plane_y(ray({0, 5, 0}, {1, 6, 0}), 0).has_value());
}

// ---- objectives ----

TEST(objectives_read_the_world) {
    World w;
    const int a = w.add(Body::solid_sphere(0.2f, 1, {1, 1, 1})), b = w.add(Body::solid_sphere(0.2f, 1, {1, 3, 1}));
    const int both[] = {a, b};
    CHECK(count_below(w, both, 2.0f) == 1);
    CHECK(inside(w, a, {0, 0, 0}, {2, 2, 2}));
    CHECK(!inside(w, a, {0, 0, 1.1f}, {2, 2, 2}));  // z counts too
    CHECK(all_slow(w, both));
    w.bodies[1].w = {0, 0, 1};
    CHECK(!all_slow(w, both));
    CHECK_NEAR(top_of(w, both), 3.2, 1e-5);
    CHECK_NEAR(top_of(w, std::span<const int>{}), kGroundY, 1e-6);
    Hold h;
    for (int i = 0; i < 59; ++i) CHECK(!h.update(true, kTimeStep, 0.5));
    CHECK(!h.update(false, kTimeStep, 0.5));
    CHECK(h.held == 0);
}

// ---- the life cycle ----

TEST(registry_lists_the_five_levels_in_menu_order) {
    CHECK(level_count() == 5);
    const char* names[] = {"Knock It Down", "Bridge Builder", "Pendulum Smash", "Chain Reaction", "Impossible Tower"};
    for (int i = 0; i < level_count(); ++i) {
        CHECK(std::strcmp(level_entry(i).name, names[i]) == 0);
        const std::unique_ptr<Level> level = make_level(i);
        CHECK(level != nullptr && std::strcmp(level->name(), names[i]) == 0);
        CHECK(level->phase == (level->has_setup() ? Phase::Setup : Phase::Running));
        CHECK(level->world.bodies[0].type == BodyType::Static);  // the ground, always first
        CHECK(!level->help().empty());
    }
    CHECK(make_level(-1) == nullptr && make_level(5) == nullptr);
}

TEST(setup_is_frozen_and_stop_returns_to_the_same_setup) {
    Pendulum p;
    p.set_pull(1.0f, 0.2f);
    const Vec3 before = p.world.bodies[static_cast<size_t>(p.bob)].pos;
    run(p, 240);
    CHECK(p.sim_time == 0);
    CHECK(p.world.bodies[static_cast<size_t>(p.bob)].pos.y == before.y);
    p.start();
    run(p, 120);
    CHECK_NEAR(p.sim_time, 1.0, 1e-6);
    CHECK(distance(p.world.bodies[static_cast<size_t>(p.bob)].pos, before) > 1.0f);
    p.stop();
    CHECK(p.phase == Phase::Setup && p.sim_time == 0 && p.status == Status::Playing);
    CHECK(p.tilt == 1.0f && p.heading == 0.2f);
    const Vec3 after = p.world.bodies[static_cast<size_t>(p.bob)].pos;
    CHECK(after.x == before.x && after.y == before.y && after.z == before.z);
    p.reset();
    CHECK(p.tilt == Pendulum::kStartTilt && p.heading == 0);
}

TEST(attempts_count_runs_not_resets) {
    Pendulum p;
    p.reset();
    CHECK(p.attempt == 1);
    p.start();
    p.stop();
    CHECK(p.attempt == 2);
    p.reset();
    CHECK(p.attempt == 2);
    p.start();
    p.reset();
    CHECK(p.attempt == 3);
}

TEST(setup_cannot_be_edited_during_a_run) {
    Pendulum p;
    p.start();
    p.set_pull(1.4f, 0.5f);
    CHECK(p.tilt == Pendulum::kStartTilt);
    Tower t;
    CHECK(t.place(3, 0, 0));
    t.start();
    run(t, 30);
    CHECK(!t.place(4, 0, 0));
    const Vec3 at = t.world.bodies[static_cast<size_t>(t.pieces[3].body)].pos;
    t.press(ray(at + Vec3{0, 5, 5}, at));  // grabbing mid-run does nothing
    CHECK(!t.holding());
    CHECK_NEAR(t.sim_time, 0.25, 1e-6);
}

// ---- placing pieces ----

TEST(a_placed_piece_drops_until_it_rests) {
    Tower t;
    CHECK(t.place(3, 0.3f, -0.2f));
    const Body& block = t.world.bodies[static_cast<size_t>(t.pieces[3].body)];
    CHECK(block.pos.x == 0.3f && block.pos.z == -0.2f);
    CHECK(block.pos.y > kGroundY + 0.6f && block.pos.y < kGroundY + 0.6f + 0.006f);
    CHECK(t.place(4, 0.3f, -0.2f));  // lands on top of the first
    CHECK_NEAR(t.world.bodies[static_cast<size_t>(t.pieces[4].body)].pos.y, kGroundY + 1.8f, 0.012);
    CHECK_NEAR(t.height(), 2.4, 0.015);
    CHECK(t.place(0, 0, 0, 0, true));  // a plank stood on end, on top
    CHECK_NEAR(t.height(), 4.8, 0.02);
}

TEST(a_piece_is_refused_off_the_pad) {
    Tower t;
    CHECK(!t.place(3, 5, 0));
    CHECK(!t.place(0, 1.6f, 0));  // a plank lying across the pad's edge sticks out
    CHECK(t.place(0, 1.5f, 0, 3.14159265f / 2));  // turned to lie along the edge, it fits
    CHECK(!t.place(99, 0, 0) && !t.place(-1, 0, 0));
    CHECK(t.placed_count() == 1);
}

TEST(pieces_are_picked_and_dropped_by_ray) {
    Tower t;
    const Vec3 home = t.pieces[3].home;
    t.press(ray(home + Vec3{0, 6, 6}, home));
    CHECK(t.holding());
    t.drag(ray({0, 8, 8}, {0, 0, 0}));  // over the middle of the pad
    t.release(ray({0, 8, 8}, {0, 0, 0}));
    CHECK(!t.holding() && t.pieces[3].placed);
    CHECK(std::fabs(t.pieces[3].pos.x) < 0.01f && std::fabs(t.pieces[3].pos.z) < 0.01f);

    // A second piece pointed at the first's top, from the side, lands on it: on the spot the ray meets
    // first, not where the same ray would have reached the ground (1.4 m further back).
    const Vec3 home2 = t.pieces[4].home;
    t.press(ray(home2 + Vec3{0, 6, 6}, home2));
    t.drag(ray({0, 8, 8}, {0, 1.2f, 0}));
    t.release(ray({0, 8, 8}, {0, 1.2f, 0}));
    CHECK_NEAR(t.height(), 2.4, 0.02);
    CHECK(std::fabs(t.pieces[4].pos.z) < 0.05f);

    // Dropped off the pad: back where it was. Cancelled: the same.
    const Vec3 was = t.pieces[4].pos;
    t.press(ray(was + Vec3{0, 6, 6}, was));
    t.release(ray({6, 8, 8}, {6, 0, 6}));
    CHECK(t.pieces[4].pos.x == was.x && t.pieces[4].pos.y == was.y);
    t.press(ray(was + Vec3{0, 6, 6}, was));
    t.drag(ray({1, 8, 8}, {1, 0, 1}));
    CHECK(t.cancel() && !t.holding());
    CHECK(t.pieces[4].pos.x == was.x && t.pieces[4].pos.y == was.y);
    CHECK(!t.cancel());
    t.unplace(4);
    CHECK(!t.pieces[4].placed && t.placed_count() == 1);
}

TEST(turning_and_tipping_a_held_piece) {
    Tower t;
    t.rotate(1);
    t.tip();
    CHECK(t.pieces[0].yaw == 0 && !t.pieces[0].upright);  // nothing held
    const Vec3 home = t.pieces[0].home;
    t.press(ray(home + Vec3{0, 6, 6}, home));
    t.drag(ray({0, 8, 8}, {0, 0, 0}));
    t.rotate(1);
    t.tip();
    t.release(ray({0, 8, 8}, {0, 0, 0}));
    CHECK(t.pieces[0].placed && t.pieces[0].upright);
    CHECK_NEAR(t.pieces[0].yaw, 3.14159265 / 12, 1e-6);
    CHECK_NEAR(t.height(), 2.4, 0.01);  // the plank stands on end
}

TEST(starting_a_run_drops_a_held_piece_back) {
    Tower t;
    t.press(ray(t.pieces[3].home + Vec3{0, 6, 6}, t.pieces[3].home));
    t.drag(ray({0, 8, 8}, {0, 0, 0}));
    t.start();
    CHECK(!t.holding() && t.placed_count() == 0);
}
