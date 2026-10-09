// The level life cycle and the piece-placing rules, exercised through real levels.
#include "test.hpp"
#include <levels/levels.hpp>
#include <levels/pendulum.hpp>
#include <levels/tower.hpp>

#include <cstring>

using namespace phys;
using namespace puzzle;

namespace {
void run(Level& l, int steps) {
    for (int i = 0; i < steps; ++i) l.step();
}
}  // namespace

TEST(registry_lists_the_five_levels_in_menu_order) {
    CHECK(level_count() == 5);
    const char* names[] = {"Knock It Down", "Bridge Builder", "Pendulum Smash", "Chain Reaction", "Impossible Tower"};
    for (int i = 0; i < level_count(); ++i) {
        CHECK(std::strcmp(level_entry(i).name, names[i]) == 0);
        const std::unique_ptr<Level> level = make_level(i);
        CHECK(level != nullptr);
        CHECK(std::strcmp(level->name(), names[i]) == 0);
        CHECK(level->status == Status::Playing);
        CHECK(level->sim_time == 0);
        CHECK(!level->help().empty());
        CHECK(level->phase == (level->has_setup() ? Phase::Setup : Phase::Running));
    }
    CHECK(make_level(-1) == nullptr);
    CHECK(make_level(5) == nullptr);
}

TEST(setup_phase_is_frozen) {
    Pendulum p;
    const Vec2 before = p.world.bodies[static_cast<size_t>(p.bob)].pos;
    run(p, 240);
    CHECK(p.phase == Phase::Setup);
    CHECK(p.sim_time == 0);
    CHECK(p.world.bodies[static_cast<size_t>(p.bob)].pos.x == before.x);
    CHECK(p.world.bodies[static_cast<size_t>(p.bob)].pos.y == before.y);
}

TEST(start_runs_and_stop_returns_to_the_same_setup) {
    Pendulum p;
    p.set_angle(-1.0f);
    const Vec2 before = p.world.bodies[static_cast<size_t>(p.bob)].pos;
    p.start();
    CHECK(p.phase == Phase::Running);
    run(p, 120);
    CHECK_NEAR(p.sim_time, 1.0, 1e-6);
    CHECK(distance(p.world.bodies[static_cast<size_t>(p.bob)].pos, before) > 1.0f);  // it swung
    p.stop();
    CHECK(p.phase == Phase::Setup);
    CHECK(p.sim_time == 0);
    CHECK(p.status == Status::Playing);
    CHECK(p.angle == -1.0f);  // the player's choice survives
    CHECK(p.world.bodies[static_cast<size_t>(p.bob)].pos.x == before.x);
    CHECK(p.world.bodies[static_cast<size_t>(p.bob)].pos.y == before.y);
    CHECK(p.world.bodies[static_cast<size_t>(p.bob)].vel.length() == 0);
    CHECK(p.world.contacts().empty());
}

TEST(stop_clears_a_finished_run) {
    Pendulum p;
    p.start();
    run(p, 120 * 10);  // the default angle is too gentle
    CHECK(p.status == Status::Failed);
    p.stop();
    CHECK(p.status == Status::Playing);
    CHECK(p.phase == Phase::Setup);
}

TEST(reset_forgets_the_setup_too) {
    Pendulum p;
    p.set_angle(-1.2f);
    p.start();
    run(p, 60);
    p.reset();
    CHECK(p.phase == Phase::Setup);
    CHECK(p.angle == Pendulum::kStartAngle);
    CHECK(p.sim_time == 0);
    const Pendulum fresh;
    CHECK(p.world.bodies.size() == fresh.world.bodies.size());
    for (size_t i = 0; i < p.world.bodies.size(); ++i) {
        CHECK(p.world.bodies[i].pos.x == fresh.world.bodies[i].pos.x);
        CHECK(p.world.bodies[i].pos.y == fresh.world.bodies[i].pos.y);
    }
}

TEST(attempts_count_runs_not_resets) {
    Pendulum p;
    CHECK(p.attempt == 1);
    p.reset();
    p.stop();
    CHECK(p.attempt == 1);  // nothing was tried yet
    p.start();
    CHECK(p.attempt == 1);
    p.stop();
    CHECK(p.attempt == 2);  // the next run will be the second
    p.reset();
    CHECK(p.attempt == 2);
    p.start();
    p.reset();
    CHECK(p.attempt == 3);
}

TEST(setup_cannot_be_edited_during_a_run) {
    Pendulum p;
    p.start();
    p.set_angle(-1.5f);
    CHECK(p.angle == Pendulum::kStartAngle);
    p.press(p.world.bodies[static_cast<size_t>(p.bob)].pos);
    p.drag({1.0f, 6.0f});
    p.release({1.0f, 6.0f});
    CHECK(p.angle == Pendulum::kStartAngle);

    Tower t;
    CHECK(t.place(3, {10, 3}));
    t.start();
    CHECK(!t.place(4, {10, 5}));
    CHECK(t.placed_count() == 1);
    run(t, 30);
    // Grabbing a piece mid-run must not pick it up (that would rebuild the world under the simulation).
    const Vec2 at = t.world.bodies[static_cast<size_t>(t.pieces[3].body)].pos;
    t.press(at);
    CHECK(!t.holding());
    t.drag({14, 6});
    t.release({14, 6});
    CHECK(t.phase == Phase::Running);
    CHECK_NEAR(t.sim_time, 0.25, 1e-6);
    CHECK(distance(t.world.bodies[static_cast<size_t>(t.pieces[3].body)].pos, at) < 0.01f);
    t.press(t.pieces[4].home);  // nor one from the tray
    CHECK(!t.holding());
}

// ---- placing pieces ----

TEST(a_placed_piece_drops_until_it_rests) {
    Tower t;
    CHECK(t.place(3, {10, 5}));  // a 1.2 m block let go 4.5 m up
    const Body& block = t.world.bodies[static_cast<size_t>(t.pieces[3].body)];
    CHECK(block.pos.x == 10.0f);
    CHECK(block.pos.y > kFloorTop + 0.6f);          // not sunk into the floor
    CHECK(block.pos.y < kFloorTop + 0.6f + 0.005f);  // within a few millimetres of it
    CHECK(t.place(4, {10, 6}));                     // and the next one lands on top of it
    CHECK_NEAR(t.world.bodies[static_cast<size_t>(t.pieces[4].body)].pos.y, kFloorTop + 1.8f, 0.01);
    CHECK_NEAR(t.height(), 2.4, 0.01);
}

TEST(a_piece_is_refused_where_it_cannot_go) {
    Tower t;
    CHECK(!t.place(3, {3, 3}));  // over the tray, outside the build area
    CHECK(t.placed_count() == 0);
    CHECK(t.place(3, {10, 3}));
    const size_t bodies = t.world.bodies.size();
    CHECK(!t.place(4, {10.5f, 1.2f}));  // inside the first block
    CHECK(t.placed_count() == 1);
    CHECK(t.world.bodies.size() == bodies);
    CHECK(!t.pieces[4].placed);
    CHECK(!t.place(99, {10, 3}));
    CHECK(!t.place(-1, {10, 3}));
}

TEST(pieces_can_be_dragged_from_the_tray_and_put_back) {
    Tower t;
    const Vec2 home = t.pieces[3].home;
    t.press(home);
    CHECK(t.holding());
    t.drag({10, 4});
    t.release({10, 4});
    CHECK(!t.holding());
    CHECK(t.pieces[3].placed);
    CHECK_NEAR(t.pieces[3].pos.y, kFloorTop + 0.6f, 0.01);

    // Let go somewhere invalid: it goes back to where it was picked up.
    const Vec2 placed_at = t.pieces[3].pos;
    t.press(placed_at);
    t.drag({3, 3});
    t.release({3, 3});
    CHECK(t.pieces[3].placed);
    CHECK(t.pieces[3].pos.x == placed_at.x && t.pieces[3].pos.y == placed_at.y);

    // Cancel mid-drag does the same.
    t.press(placed_at);
    t.drag({14, 6});
    CHECK(t.cancel());
    CHECK(!t.holding());
    CHECK(t.pieces[3].pos.x == placed_at.x && t.pieces[3].pos.y == placed_at.y);
    CHECK(!t.cancel());

    t.unplace(3);
    CHECK(!t.pieces[3].placed);
    CHECK(t.placed_count() == 0);
    CHECK(t.world.bodies.size() == Tower().world.bodies.size());
}

TEST(a_tray_piece_let_go_over_the_tray_stays_in_the_tray) {
    Tower t;
    t.press(t.pieces[0].home);
    t.release(t.pieces[0].home);
    CHECK(!t.pieces[0].placed);
    CHECK(t.placed_count() == 0);
}

TEST(rotating_a_held_piece) {
    Tower t;
    t.rotate(1);  // nothing held: nothing happens
    CHECK(t.pieces[0].angle == 0);
    t.press(t.pieces[0].home);
    t.drag({10, 5});
    t.rotate(1);
    t.release({10, 5});
    CHECK(t.pieces[0].placed);
    CHECK_NEAR(t.pieces[0].angle, kPi / 2, 1e-5);
    CHECK_NEAR(t.height(), 2.4, 0.01);  // the plank stands on end
}

TEST(starting_a_run_drops_a_held_piece_back) {
    Tower t;
    t.press(t.pieces[3].home);
    t.drag({10, 5});
    t.start();
    CHECK(!t.holding());
    CHECK(t.placed_count() == 0);
}
