// Each level played headless: arrange, start, step until the level decides.
#include "test.hpp"
#include <levels/bridge.hpp>
#include <levels/chain_reaction.hpp>
#include <levels/knock_down.hpp>
#include <levels/pendulum.hpp>
#include <levels/tower.hpp>

#include <cmath>

using namespace phys3d;
using namespace puzzle3d;

namespace {

constexpr Real kDeg = 3.14159265f / 180.0f;

Status play(Level& level, double max_seconds = 40.0) {
    level.start();
    const int steps = static_cast<int>(max_seconds / static_cast<double>(kTimeStep));
    for (int i = 0; i < steps && level.status == Status::Playing; ++i) level.step();
    return level.status;
}
void run(Level& l, double seconds) {
    for (int i = 0; i < static_cast<int>(std::lround(seconds / static_cast<double>(kTimeStep))); ++i) l.step();
}
// Fires one ball at each aim, waiting 1.5 s between them, then plays on until the level decides.
void fire_all(KnockDown& k, std::initializer_list<Vec3> aims) {
    auto next = aims.begin();
    double wait_until = 0;
    for (int i = 0; i < 120 * 30 && k.status == Status::Playing; ++i) {
        if (next != aims.end() && k.sim_time >= wait_until && k.fire(*next)) {
            ++next;
            wait_until = k.sim_time + 1.5;
        }
        k.step();
    }
}
bool same_bodies(const World& a, const World& b) {
    if (a.bodies.size() != b.bodies.size() || a.joints.size() != b.joints.size()) return false;
    for (size_t i = 0; i < a.bodies.size(); ++i) {
        const Body &p = a.bodies[i], &q = b.bodies[i];
        if (p.pos.x != q.pos.x || p.pos.y != q.pos.y || p.pos.z != q.pos.z || p.q.w != q.q.w || p.vel.length() != q.vel.length()) return false;
    }
    return true;
}
constexpr Vec3 kWin1{12, 4.5f, 0.4f}, kWin2{11, 4.5f, -0.4f}, kWin3{11, 4.5f, 0};  // found by searching

}  // namespace

// ---- Level 1: Knock It Down ----

TEST(knock_down_pyramid_stands_and_balls_are_real_bodies) {
    KnockDown k;
    CHECK(k.targets.size() == 14);
    const KnockDown start;
    run(k, 5.0);
    CHECK(k.fallen == 0 && k.status == Status::Playing);
    for (int t : k.targets) CHECK(distance(k.world.bodies[static_cast<size_t>(t)].pos, start.world.bodies[static_cast<size_t>(t)].pos) < 0.02f);
    const size_t before = k.world.bodies.size();
    CHECK(k.fire(kWin1));
    CHECK(k.world.bodies.size() == before + 1);
    const Body& ball = k.world.bodies.back();
    CHECK(ball.type == BodyType::Dynamic);
    const Vec3 v = k.launch_velocity(kWin1);
    CHECK(ball.vel.x == v.x && ball.vel.y == v.y && ball.vel.z == v.z);
}

TEST(knock_down_launch_speed_is_clamped_and_points_at_the_aim) {
    const KnockDown k;
    CHECK_NEAR(k.launch_velocity(KnockDown::kMuzzle + Vec3{0.1f, 0, 0}).length(), KnockDown::kMinSpeed, 1e-4);
    const Vec3 far = k.launch_velocity(KnockDown::kMuzzle + Vec3{30, 20, 10});
    CHECK_NEAR(far.length(), KnockDown::kMaxSpeed, 1e-3);
    CHECK_NEAR(far.y / far.x, 20.0 / 30.0, 1e-4);
    CHECK_NEAR(far.z / far.x, 10.0 / 30.0, 1e-4);
    CHECK_NEAR(k.launch_velocity(KnockDown::kMuzzle + Vec3{0, 0, 5}).z, 5 * KnockDown::kSpeedPerMeter, 1e-4);
    // Aiming by ray: at the crates it is the point where the ray meets the first crate.
    const Vec3 crate = k.world.bodies[static_cast<size_t>(k.targets.back())].pos;  // the top one
    const Vec3 from{0, crate.y, 0};
    const Vec3 aim = k.aim_point({from, (crate - from).normalized()});
    CHECK(std::fabs(aim.x - (crate.x - 0.35f)) < 1e-3f);
}

TEST(knock_down_shot_limit_and_reload) {
    KnockDown k;
    CHECK(k.fire(kWin1));
    CHECK(!k.fire(kWin1));  // a double click is one shot
    run(k, 0.5);
    for (int i = 1; i < KnockDown::kShots; ++i) {
        CHECK(k.fire({0, 2, -8}));
        run(k, 0.5);
    }
    CHECK(k.shots_left == 0 && !k.can_fire() && !k.fire(kWin1));
}

TEST(knock_down_three_well_placed_balls_win_and_wasted_ones_fail) {
    KnockDown k;
    fire_all(k, {kWin1, kWin2, kWin3});
    CHECK(k.status == Status::Won);
    CHECK(k.shots_used() == 3);
    for (int t : k.targets) CHECK(k.world.bodies[static_cast<size_t>(t)].pos.y < KnockDown::kFallLine);

    KnockDown same;  // the same aim four times leaves some standing
    fire_all(same, {kWin1, kWin1, kWin1, kWin1});
    CHECK(same.status == Status::Failed);
    CHECK(same.fallen < 14);

    KnockDown wasted;
    fire_all(wasted, {{0, 2, -8}, {0, 2, -8}, {0, 2, -8}, {0, 2, -8}});
    CHECK(wasted.status == Status::Failed && wasted.fallen == 0);
}

TEST(knock_down_reset_and_replay) {
    KnockDown a;
    fire_all(a, {kWin1, kWin2, kWin3});
    a.reset();
    const KnockDown fresh;
    CHECK(same_bodies(a.world, fresh.world));
    CHECK(a.shots_left == KnockDown::kShots && a.sim_time == 0 && a.status == Status::Playing && a.fallen == 0);
    fire_all(a, {kWin1, kWin2, kWin3});
    KnockDown b;
    fire_all(b, {kWin1, kWin2, kWin3});
    CHECK(a.win_time == b.win_time && same_bodies(a.world, b.world));
}

// ---- Level 3: Pendulum Smash ----

TEST(pendulum_swings_on_a_real_rod) {
    Pendulum p;
    CHECK(p.world.joints.size() == 1 && p.world.joints[0].type == JointType::Distance);
    p.set_pull(60 * kDeg, 0);
    p.start();
    for (int i = 0; i < 60; ++i) {
        p.step();
        CHECK_NEAR(distance(p.world.bodies[static_cast<size_t>(p.bob)].pos, Pendulum::kPivot), Pendulum::kRod, 0.01);
    }
}

TEST(pendulum_needs_the_right_pull_and_the_right_aim) {
    Pendulum gentle;  // the starting pull
    CHECK(play(gentle) == Status::Failed);
    for (const Real tilt : {50.0f, 56.0f, 62.0f}) {
        Pendulum p;
        p.set_pull(tilt * kDeg, 0);
        CHECK(play(p) == Status::Won);
        CHECK(p.target_in_bin());
    }
    Pendulum hard;
    hard.set_pull(80 * kDeg, 0);
    CHECK(play(hard) == Status::Failed);
    CHECK(hard.world.bodies[static_cast<size_t>(hard.target)].pos.x > Pendulum::kZoneHi.x);  // over the bin
    for (const Real off : {-10.0f, 10.0f}) {  // a good pull, aimed ten degrees wide
        Pendulum p;
        p.set_pull(60 * kDeg, off * kDeg);
        CHECK(play(p) == Status::Failed);
        CHECK(!p.target_in_bin());
    }
}

TEST(pendulum_is_pulled_by_dragging_on_the_sphere) {
    Pendulum p;
    const Vec3 bob = p.world.bodies[static_cast<size_t>(p.bob)].pos;
    const Vec3 eye{-12, 6, 10};
    p.press({eye, (bob - eye).normalized()});
    const Vec3 to = Pendulum::bob_position(1.0f, 0.3f);
    p.drag({eye, (to - eye).normalized()});
    CHECK_NEAR(p.tilt, 1.0, 1e-3);
    CHECK_NEAR(p.heading, 0.3, 1e-3);
    CHECK(distance(p.world.bodies[static_cast<size_t>(p.bob)].pos, to) < 0.01f);
    CHECK(p.cancel());
    CHECK(p.tilt == Pendulum::kStartTilt && p.heading == 0);
    p.press({eye, (Vec3{12, 1, 6} - eye).normalized()});  // nowhere near the weight
    p.drag({eye, (to - eye).normalized()});
    CHECK(p.tilt == Pendulum::kStartTilt);
    p.rotate(1);
    CHECK_NEAR(p.heading, kDeg, 1e-6);
    p.set_pull(3.0f, 2.0f);
    CHECK(p.tilt == Pendulum::kMaxTilt);
}

// ---- Level 4: Chain Reaction ----

namespace {
void fill_gap(ChainReaction& c, Real yaw_error = 0) {
    for (int i = 0; i < ChainReaction::kGapCount; ++i) {
        const int k = ChainReaction::kGapFirst + i;
        const Vec3 p = ChainReaction::slot_position(k);
        CHECK(c.place(i, p.x, p.z, ChainReaction::slot_yaw(k) + yaw_error));
    }
}
}  // namespace

TEST(chain_with_the_gap_left_open_stops) {
    ChainReaction c;
    const Vec3 target = c.world.bodies[static_cast<size_t>(c.target)].pos;
    CHECK(play(c) == Status::Failed);
    CHECK(distance(c.world.bodies[static_cast<size_t>(c.target)].pos, target) < 0.01f);  // never reached
}

TEST(chain_spares_turned_to_follow_the_curve_win) {
    ChainReaction c;
    fill_gap(c);
    CHECK(play(c) == Status::Won);
    CHECK(c.target_in_bin() && c.win_time > 3.0);
    ChainReaction rough;  // turned only to the nearest 15 degrees, as Q and E do
    for (int i = 0; i < ChainReaction::kGapCount; ++i) {
        const int k = ChainReaction::kGapFirst + i;
        const Vec3 p = ChainReaction::slot_position(k);
        const Real step = 3.14159265f / 12;
        CHECK(rough.place(i, p.x, p.z, std::round(ChainReaction::slot_yaw(k) / step) * step));
    }
    CHECK(play(rough) == Status::Won);
}

TEST(chain_spares_left_unturned_or_too_few_fail) {
    ChainReaction unturned;
    for (int i = 0; i < ChainReaction::kGapCount; ++i) {
        const Vec3 p = ChainReaction::slot_position(ChainReaction::kGapFirst + i);
        CHECK(unturned.place(i, p.x, p.z, 0));
    }
    CHECK(play(unturned) == Status::Failed);
    ChainReaction askew;
    fill_gap(askew, 60 * kDeg);
    CHECK(play(askew) == Status::Failed);
    ChainReaction one;
    const int k = ChainReaction::kGapFirst + 1;
    CHECK(one.place(0, ChainReaction::slot_position(k).x, ChainReaction::slot_position(k).z, ChainReaction::slot_yaw(k)));
    CHECK(play(one) == Status::Failed);
}

TEST(chain_spares_only_go_in_the_gap_and_not_on_anything) {
    ChainReaction c;
    CHECK(!c.place(0, 0, 0));  // the middle of the circle: outside the gap
    const Vec3 neighbour = ChainReaction::slot_position(ChainReaction::kGapFirst - 1);
    CHECK(!c.place(0, neighbour.x, neighbour.z));  // on top of the domino beside the gap
    fill_gap(c);
    const Vec3 p = ChainReaction::slot_position(ChainReaction::kGapFirst);
    CHECK_NEAR(c.pieces[0].pos.y, kGroundY + ChainReaction::kDominoHalf.y, 0.006);
    c.unplace(0);
    CHECK(c.placed_count() == 2);
    const Vec3 q = ChainReaction::slot_position(ChainReaction::kGapFirst + 1);
    CHECK(!c.place(0, q.x, q.z));  // on top of another spare: too high for the gap
    CHECK(c.placed_count() == 2);
    CHECK(c.place(0, p.x, p.z, ChainReaction::slot_yaw(ChainReaction::kGapFirst)));
}

TEST(chain_reset_and_replay) {
    ChainReaction a;
    fill_gap(a);
    CHECK(play(a) == Status::Won);
    const double t = a.win_time;
    a.reset();
    CHECK(a.placed_count() == 0 && same_bodies(a.world, ChainReaction().world));
    fill_gap(a);
    CHECK(play(a) == Status::Won && a.win_time == t);
}

// ---- Level 5: Impossible Tower ----

TEST(tower_nothing_or_too_short_fails) {
    Tower empty;
    CHECK(play(empty) == Status::Failed);
    Tower low;
    CHECK(low.place(3, 0, 0) && low.place(4, 0, 0));
    CHECK(play(low) == Status::Failed && low.best_height == 0);
}

TEST(tower_wins_only_after_standing_for_the_hold_time) {
    Tower t;
    CHECK(t.place(3, 0, 0) && t.place(0, 0, 0, 0, true) && t.place(4, 0, 0));
    CHECK(t.height() > Tower::kGoal);
    t.start();
    while (t.sim_time < Tower::kHoldTime - 0.05) {
        t.step();
        CHECK(t.status == Status::Playing);
    }
    run(t, 2.0);
    CHECK(t.status == Status::Won && t.win_time >= Tower::kHoldTime);
    CHECK_NEAR(t.best_height, t.height(), 1e-4);
}

TEST(tower_tall_at_the_start_but_toppling_fails) {
    Tower t;
    CHECK(t.place(3, -1, 0) && t.place(4, -0.4f, 0) && t.place(0, 0.2f, 0) && t.place(1, 0.8f, 0));
    CHECK(t.place(2, 1.4f, 0, 3.14159265f / 2) && t.place(6, 1.5f, 0));
    CHECK(t.height() > Tower::kGoal);
    CHECK(play(t) == Status::Failed);
    CHECK(t.height() < Tower::kGoal && t.best_height == 0);
}

TEST(tower_best_height_is_kept_and_only_grows) {
    Tower t;
    for (const int i : {3, 4, 0, 1, 2, 6, 5}) CHECK(t.place(i, 0, 0));  // everything flat: about 5.4 m
    CHECK(play(t) == Status::Won);
    const Real best = t.best_height;
    CHECK(best > 5.3f);
    t.reset();
    CHECK(t.best_height == best && t.placed_count() == 0);
    CHECK(t.place(0, 0, 0, 0, true) && t.place(1, 0, 0, 0, true));  // 4.8 m: passes, but lower
    CHECK(play(t) == Status::Won && t.best_height == best);
}

// ---- Level 2: Bridge Builder ----

namespace {
struct Edge {
    int c1, r1, c2, r2;
};
void add(Bridge& b, std::initializer_list<Edge> edges) {
    for (const Edge& e : edges) CHECK(b.add_beam(Bridge::node(e.c1, e.r1), Bridge::node(e.c2, e.r2)));
}
void deck_and_ends(Bridge& b) {
    add(b, {{0, 1, 1, 1}, {1, 1, 2, 1}, {2, 1, 3, 1}, {3, 1, 4, 1}, {0, 0, 1, 1}, {4, 0, 4, 1}, {4, 0, 3, 1}});
}
}  // namespace

TEST(bridge_with_no_bridge_or_a_bare_deck_the_car_falls) {
    Bridge none;
    CHECK(play(none) == Status::Failed && none.chassis().pos.y < Bridge::kFallY);
    Bridge deck;
    add(deck, {{0, 1, 1, 1}, {1, 1, 2, 1}, {2, 1, 3, 1}, {3, 1, 4, 1}});
    CHECK(play(deck) == Status::Failed && deck.chassis().pos.y < Bridge::kFallY);
}

TEST(bridge_seven_or_nine_beams_are_not_enough) {
    Bridge seven;
    deck_and_ends(seven);
    CHECK(play(seven) == Status::Failed && seven.chassis().pos.x < Bridge::kFinishX);
    Bridge nine;  // a left triangle, but no post on the right
    add(nine, {{0, 1, 1, 1}, {1, 1, 2, 1}, {2, 1, 3, 1}, {3, 1, 4, 1}, {0, 0, 1, 1}, {4, 0, 3, 1}, {0, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 2, 1}});
    CHECK(play(nine) == Status::Failed);
}

TEST(bridge_triangulated_designs_carry_the_car) {
    Bridge left;
    deck_and_ends(left);
    add(left, {{0, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 2, 1}});
    CHECK(play(left) == Status::Won);
    CHECK(left.chassis().pos.x >= Bridge::kFinishX && left.win_time > 2.0);
    run(left, 2.0);
    CHECK(left.chassis().pos.x < Bridge::kFinishX + 1.5f);  // it braked
    Bridge centre;  // another answer: a V under the middle
    deck_and_ends(centre);
    add(centre, {{1, 1, 2, 0}, {2, 0, 3, 1}, {2, 0, 2, 1}});
    CHECK(play(centre) == Status::Won);
    // The hinges held: no joint came apart by more than a millimetre.
    for (const Joint& j : centre.world.joints) CHECK(distance(j.world_anchor_a(centre.world.bodies), j.world_anchor_b(centre.world.bodies)) < 0.001f);
}

TEST(bridge_beams_are_hinged_bodies_and_follow_the_rules) {
    Bridge b;
    const size_t bodies = b.world.bodies.size(), joints = b.world.joints.size();
    CHECK(b.add_beam(Bridge::node(0, 1), Bridge::node(1, 1)));
    CHECK(b.world.bodies.size() == bodies + 1 && b.world.joints.size() == joints + 1);
    CHECK(b.world.joints.back().type == JointType::Hinge && b.world.joints.back().a == -1);
    CHECK(b.add_beam(Bridge::node(1, 1), Bridge::node(2, 1)));
    CHECK(b.world.joints.back().a >= 0 && b.world.joints.back().b >= 0);
    CHECK(!b.add_beam(Bridge::node(2, 1), Bridge::node(1, 1)));  // already there
    CHECK(!b.add_beam(Bridge::node(1, 1), Bridge::node(3, 1)));  // not neighbours
    CHECK(!b.add_beam(Bridge::node(0, 0), Bridge::node(0, 1)));  // rock to rock
    CHECK(!Bridge::is_anchor(Bridge::node(4, 1)) && Bridge::is_anchor(Bridge::node(4, 0)));
    CHECK(b.remove_beam(Bridge::node(1, 1), Bridge::node(2, 1)) && b.beams.size() == 1);
    int added = static_cast<int>(b.beams.size());
    for (int c = 0; c < Bridge::kColumns - 1; ++c)
        for (int r = 0; r < Bridge::kRows; ++r) {
            added += b.add_beam(Bridge::node(c, r), Bridge::node(c + 1, r)) ? 1 : 0;
            added += b.add_beam(Bridge::node(c, r), Bridge::node(c, (r + 1) % 2)) ? 1 : 0;
            added += b.add_beam(Bridge::node(c, r), Bridge::node(c + 1, (r + 1) % 2)) ? 1 : 0;
        }
    CHECK(added == Bridge::kMaxBeams && static_cast<int>(b.beams.size()) == Bridge::kMaxBeams);
}

TEST(bridge_is_built_by_ray_and_frozen_during_a_run) {
    Bridge b;
    const Vec3 eye{0, 4, 16};
    auto at = [&](Vec3 p) { return Ray{eye, (p - eye).normalized()}; };
    const Vec3 p = Bridge::node_position(Bridge::node(1, 1)), q = Bridge::node_position(Bridge::node(2, 1));
    b.press(at(p + Vec3{0.2f, 0.1f, 0}));
    b.release(at(q));
    CHECK(b.beams.size() == 1);
    b.press(at(p));
    b.release(at({0, 9, 0}));  // let go over nothing
    CHECK(b.beams.size() == 1);
    b.press(at(p));
    CHECK(b.cancel());
    b.release(at(Bridge::node_position(Bridge::node(1, 0))));
    CHECK(b.beams.size() == 1);
    CHECK(!b.cancel());
    b.start();
    b.press(at((p + q) * 0.5f));
    b.release(at((p + q) * 0.5f));
    CHECK(b.beams.size() == 1);
    CHECK(!b.add_beam(Bridge::node(2, 1), Bridge::node(3, 1)));
    CHECK(!b.remove_beam(Bridge::node(1, 1), Bridge::node(2, 1)));
    CHECK(b.beams.size() == 1);
    b.stop();
    b.press(at((p + q) * 0.5f));  // a click on a beam takes it away
    b.release(at((p + q) * 0.5f));
    CHECK(b.beams.empty());
}

TEST(bridge_reset_and_replay) {
    Bridge a;
    deck_and_ends(a);
    add(a, {{0, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 2, 1}});
    CHECK(play(a) == Status::Won);
    const double t = a.win_time;
    a.stop();
    CHECK(a.beams.size() == 10 && a.chassis().pos.x < Bridge::kLeftEdge);
    CHECK(play(a) == Status::Won && a.win_time == t);
    a.reset();
    CHECK(a.beams.empty() && same_bodies(a.world, Bridge().world));
}
