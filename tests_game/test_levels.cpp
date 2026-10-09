#include "test.hpp"
#include <levels/knock_down.hpp>

#include <cmath>

using namespace phys;
using namespace puzzle;

namespace {

constexpr Vec2 kGoodAim{14.5f, 6.0f};  // an arc through the middle of the pyramid at full speed
constexpr Vec2 kBadAim{0.6f, 2.6f};    // backwards into the wall behind the launcher

void run(KnockDown& k, double seconds) {
    const int steps = static_cast<int>(std::lround(seconds / static_cast<double>(kTimeStep)));
    for (int i = 0; i < steps; ++i) k.step();
}
// Fires every shot at `aim` as soon as the launcher allows, then plays on until the level decides.
void play_out(KnockDown& k, Vec2 aim) {
    for (int i = 0; i < 120 * 30 && k.status == Status::Playing; ++i) {
        k.fire(aim);
        k.step();
    }
}
bool same_state(const KnockDown& a, const KnockDown& b) {
    if (a.world.bodies.size() != b.world.bodies.size() || a.world.joints.size() != b.world.joints.size()) return false;
    for (size_t i = 0; i < a.world.bodies.size(); ++i) {
        const Body &p = a.world.bodies[i], &q = b.world.bodies[i];
        if (p.pos.x != q.pos.x || p.pos.y != q.pos.y || p.angle != q.angle) return false;
        if (p.vel.x != q.vel.x || p.vel.y != q.vel.y || p.w != q.w || p.awake != q.awake) return false;
    }
    return a.shots_left == b.shots_left && a.fallen == b.fallen && a.sim_time == b.sim_time &&
           a.last_shot_time == b.last_shot_time && a.win_time == b.win_time && a.status == b.status &&
           a.targets == b.targets && a.first_ball == b.first_ball;
}

}  // namespace

TEST(knock_down_stack_stands_on_its_own) {
    KnockDown k;
    const KnockDown start;
    run(k, 5.0);
    CHECK(k.status == Status::Playing);
    CHECK(k.fallen == 0);
    for (int t : k.targets)
        CHECK(distance(k.world.bodies[static_cast<size_t>(t)].pos, start.world.bodies[static_cast<size_t>(t)].pos) < 0.02f);
}

TEST(knock_down_launch_speed_is_clamped_and_points_at_the_aim) {
    const KnockDown k;
    const Vec2 near = k.launch_velocity(KnockDown::kMuzzle + Vec2{0.1f, 0});
    CHECK_NEAR(near.length(), KnockDown::kMinSpeed, 1e-4);
    const Vec2 far = k.launch_velocity(KnockDown::kMuzzle + Vec2{30.0f, 40.0f});
    CHECK_NEAR(far.length(), KnockDown::kMaxSpeed, 1e-3);
    CHECK_NEAR(far.y / far.x, 40.0 / 30.0, 1e-4);
    const Vec2 mid = k.launch_velocity(KnockDown::kMuzzle + Vec2{0, 5.0f});
    CHECK_NEAR(mid.x, 0, 1e-5);
    CHECK_NEAR(mid.y, 5.0f * KnockDown::kSpeedPerMeter, 1e-4);
    CHECK(std::isfinite(k.launch_velocity(KnockDown::kMuzzle).x));  // aiming exactly at the muzzle
}

TEST(knock_down_a_shot_is_a_real_body_with_the_launch_velocity) {
    KnockDown k;
    const size_t before = k.world.bodies.size();
    CHECK(k.fire(kGoodAim));
    CHECK(k.world.bodies.size() == before + 1);
    const Body& ball = k.world.bodies.back();
    CHECK(ball.type == BodyType::Dynamic);
    CHECK(ball.pos.x == KnockDown::kMuzzle.x && ball.pos.y == KnockDown::kMuzzle.y);
    CHECK(ball.vel.x == k.launch_velocity(kGoodAim).x && ball.vel.y == k.launch_velocity(kGoodAim).y);
    CHECK(k.shots_left == KnockDown::kShots - 1);
}

TEST(knock_down_shot_limit) {
    KnockDown k;
    const size_t before = k.world.bodies.size();
    for (int i = 0; i < KnockDown::kShots; ++i) {
        CHECK(k.fire(kBadAim));
        run(k, 0.5);
    }
    CHECK(k.shots_left == 0);
    CHECK(!k.can_fire());
    CHECK(!k.fire(kBadAim));
    CHECK(k.shots_left == 0);
    CHECK(k.world.bodies.size() == before + static_cast<size_t>(KnockDown::kShots));
}

// A held button or a double click must not empty the launcher in one frame.
TEST(knock_down_cannot_fire_again_until_reloaded) {
    KnockDown k;
    CHECK(k.fire(kGoodAim));
    CHECK(!k.fire(kGoodAim));
    run(k, static_cast<double>(KnockDown::kReloadTime) - 0.1);
    CHECK(!k.fire(kGoodAim));
    CHECK(k.shots_left == KnockDown::kShots - 1);
    run(k, 0.15);
    CHECK(k.fire(kGoodAim));
    CHECK(k.shots_left == KnockDown::kShots - 2);
}

TEST(knock_down_good_shots_win) {
    KnockDown k;
    play_out(k, kGoodAim);
    CHECK(k.status == Status::Won);
    CHECK(k.fallen == static_cast<int>(k.targets.size()));
    CHECK(k.win_time > 0 && k.win_time <= k.sim_time);
    for (int t : k.targets) CHECK(k.world.bodies[static_cast<size_t>(t)].pos.y < KnockDown::kFallLine);
    // Won is final: no more shots, and it stays won.
    CHECK(!k.fire(kGoodAim));
    const double won_at = k.win_time;
    run(k, 2.0);
    CHECK(k.status == Status::Won);
    CHECK(k.win_time == won_at);
}

TEST(knock_down_one_block_left_standing_is_not_a_win) {
    KnockDown k;
    // Put every block but one on the floor by hand; the objective must still say no.
    for (size_t i = 1; i < k.targets.size(); ++i) {
        Body& b = k.world.bodies[static_cast<size_t>(k.targets[i])];
        b.pos = {5.0f + static_cast<Real>(i), kFloorTop + 0.36f};
        b.wake();
    }
    run(k, 1.0);
    CHECK(k.fallen == static_cast<int>(k.targets.size()) - 1);
    CHECK(k.status == Status::Playing);
}

TEST(knock_down_no_firing_once_solved_even_with_balls_left) {
    KnockDown k;
    for (size_t i = 0; i < k.targets.size(); ++i) {
        Body& b = k.world.bodies[static_cast<size_t>(k.targets[i])];
        b.pos = {5.0f + static_cast<Real>(i), kFloorTop + 0.36f};
        b.wake();
    }
    k.step();
    CHECK(k.status == Status::Won);
    CHECK(k.shots_left == KnockDown::kShots);
    CHECK(!k.can_fire());
    CHECK(!k.fire(kGoodAim));
    CHECK(k.world.bodies.size() == k.first_ball);
}

TEST(knock_down_wasted_shots_fail_only_after_the_last_one) {
    KnockDown k;
    run(k, 10.0);  // waiting costs nothing while shots remain
    CHECK(k.status == Status::Playing);
    for (int i = 0; i < KnockDown::kShots - 1; ++i) {
        CHECK(k.fire(kBadAim));
        run(k, 8.0);
        CHECK(k.status == Status::Playing);
    }
    CHECK(k.fire(kBadAim));
    run(k, 0.25);
    CHECK(k.status == Status::Playing);  // the last ball gets its chance
    run(k, static_cast<double>(KnockDown::kSettleTime));
    CHECK(k.status == Status::Failed);
    CHECK(k.fallen == 0);
}

TEST(knock_down_reset_restores_everything) {
    const KnockDown fresh;
    KnockDown k;
    play_out(k, kGoodAim);
    CHECK(!same_state(k, fresh));
    k.reset();
    CHECK(same_state(k, fresh));
    CHECK(k.world.contacts().empty());
    CHECK(k.can_fire());

    KnockDown failed;
    play_out(failed, kBadAim);
    CHECK(failed.status == Status::Failed);
    failed.reset();
    CHECK(same_state(failed, fresh));
}

// Same shots at the same simulated times give the same result, whatever was played before the reset.
TEST(knock_down_replay_after_reset_is_identical) {
    KnockDown first;
    play_out(first, kGoodAim);

    KnockDown second;
    play_out(second, kBadAim);
    second.reset();
    play_out(second, kGoodAim);
    CHECK(same_state(first, second));
}

// ---------------------------------------------------------------------------------------------
// The setup levels. Each is played headless: arrange, start, step until the level decides.
#include <levels/bridge.hpp>
#include <levels/chain_reaction.hpp>
#include <levels/pendulum.hpp>
#include <levels/tower.hpp>

namespace {

constexpr Real kDeg = kPi / 180.0f;

Status play(Level& level, double max_seconds = 40.0) {
    level.start();
    const int steps = static_cast<int>(max_seconds / static_cast<double>(kTimeStep));
    for (int i = 0; i < steps && level.status == Status::Playing; ++i) level.step();
    return level.status;
}

}  // namespace

// ---- Level 3: Pendulum Smash ----

TEST(pendulum_is_a_real_jointed_weight) {
    Pendulum p;
    CHECK(p.world.joints.size() == 1);
    CHECK(p.world.joints[0].type == JointType::Distance);
    p.set_angle(-60 * kDeg);
    p.start();
    for (int i = 0; i < 90; ++i) {  // before it reaches the crate: the rod keeps its length
        p.step();
        CHECK_NEAR(distance(p.world.bodies[static_cast<size_t>(p.bob)].pos, Pendulum::kPivot), Pendulum::kRod, 0.01);
    }
    CHECK(p.world.bodies[static_cast<size_t>(p.bob)].vel.length() > 3.0f);  // gravity did that, nothing else
}

TEST(pendulum_a_gentle_swing_leaves_the_crate_on_the_table) {
    Pendulum p;  // the starting angle
    CHECK(play(p) == Status::Failed);
    CHECK(!p.crate_in_bin());
    CHECK(p.world.bodies[static_cast<size_t>(p.crate)].pos.y > 4.0f);
}

TEST(pendulum_the_right_swing_wins) {
    for (const Real degrees : {56.0f, 64.0f, 72.0f}) {
        Pendulum p;
        p.set_angle(-degrees * kDeg);
        CHECK(play(p) == Status::Won);
        CHECK(p.crate_in_bin());
        CHECK(p.win_time > 1.0);
    }
}

TEST(pendulum_too_hard_a_swing_overshoots) {
    Pendulum p;
    p.set_angle(-85 * kDeg);
    CHECK(play(p) == Status::Failed);
    CHECK(!p.crate_in_bin());
    CHECK(p.world.bodies[static_cast<size_t>(p.crate)].pos.y < 2.0f);  // it did leave the table
}

// Passing through the bin is not enough: the crate has to stay.
TEST(pendulum_win_needs_the_crate_at_rest_in_the_bin) {
    Pendulum p;
    p.start();
    Body& crate = p.world.bodies[static_cast<size_t>(p.crate)];
    crate.pos = {13.9f, 3.0f};  // dropped into the bin from above
    crate.vel = {};
    int steps = 0;
    while (p.status == Status::Playing && steps < 1200) {
        p.step();
        ++steps;
        if (p.status == Status::Won) break;
        if (p.crate_in_bin() && steps < 30) CHECK(p.status == Status::Playing);
    }
    CHECK(p.status == Status::Won);
    CHECK(p.win_time >= Pendulum::kRestTime);
}

TEST(pendulum_angle_is_clamped_and_set_by_dragging) {
    Pendulum p;
    p.set_angle(-3.0f);
    CHECK(p.angle == Pendulum::kMinAngle);
    p.set_angle(1.0f);
    CHECK(p.angle == Pendulum::kMaxAngle);

    p.set_angle(-0.5f);
    p.press(p.world.bodies[static_cast<size_t>(p.bob)].pos);
    p.drag(Pendulum::kPivot + Vec2{-3.0f, -3.0f});  // 45 degrees to the left, at any distance
    CHECK_NEAR(p.angle, -kPi / 4, 1e-5);
    const Vec2 bob = p.world.bodies[static_cast<size_t>(p.bob)].pos;
    CHECK_NEAR(distance(bob, Pendulum::kPivot), Pendulum::kRod, 1e-4);
    CHECK(p.cancel());
    CHECK_NEAR(p.angle, -0.5f, 1e-6);

    p.press({15.0f, 8.0f});  // nowhere near the weight: not a grab
    p.drag(Pendulum::kPivot + Vec2{-3.0f, -3.0f});
    CHECK_NEAR(p.angle, -0.5f, 1e-6);

    p.rotate(1);
    CHECK_NEAR(p.angle, -0.5f - kDeg, 1e-6);
}

TEST(pendulum_same_angle_same_result) {
    Pendulum a, b;
    a.set_angle(-64 * kDeg);
    b.set_angle(-80 * kDeg);
    play(b);
    b.stop();
    b.set_angle(-64 * kDeg);
    play(a);
    play(b);
    CHECK(a.win_time == b.win_time);
    CHECK(a.world.bodies[static_cast<size_t>(a.crate)].pos.x == b.world.bodies[static_cast<size_t>(b.crate)].pos.x);
}

// ---- Level 4: Chain Reaction ----

TEST(chain_without_the_missing_dominoes_stops_at_the_gap) {
    ChainReaction c;
    const Real target_x = c.world.bodies[static_cast<size_t>(c.target)].pos.x;
    CHECK(play(c) == Status::Failed);
    CHECK(!c.target_in_bin());
    CHECK_NEAR(c.world.bodies[static_cast<size_t>(c.target)].pos.x, target_x, 0.01);  // never touched
    CHECK(c.world.bodies[static_cast<size_t>(c.starter)].pos.x > 4.0f);                // but the first ball did roll
}

TEST(chain_evenly_spaced_spares_carry_it_to_the_bin) {
    ChainReaction c;
    const Real left = ChainReaction::kGapLeft, right = ChainReaction::kGapRight;
    for (int i = 0; i < 3; ++i) CHECK(c.place(i, {left + (right - left) * static_cast<Real>(i + 1) / 4.0f, 4.0f}));
    CHECK(play(c) == Status::Won);
    CHECK(c.target_in_bin());
    CHECK(c.win_time > 3.0);  // it takes that long for the chain to get there: nothing was skipped
}

TEST(chain_one_spare_is_not_enough_and_bunched_spares_do_not_reach) {
    ChainReaction one;
    CHECK(one.place(0, {8.0f, 4.0f}));
    CHECK(play(one) == Status::Failed);

    ChainReaction bunched;
    for (int i = 0; i < 3; ++i) CHECK(bunched.place(i, {ChainReaction::kGapLeft + 0.4f + 0.3f * static_cast<Real>(i), 4.0f}));
    CHECK(play(bunched) == Status::Failed);
    CHECK(!bunched.target_in_bin());
}

TEST(chain_spares_only_go_in_the_gap_and_land_on_the_shelf) {
    ChainReaction c;
    CHECK(!c.place(0, {12.6f, 4.0f}));  // beside the target ball: outside the marked gap
    CHECK(!c.place(0, {5.5f, 4.0f}));
    CHECK(c.place(0, {8.0f, 4.6f}));
    CHECK_NEAR(c.pieces[0].pos.y, ChainReaction::kShelfTop + ChainReaction::kDominoHalf.y, 0.005);
    CHECK(c.pieces[0].pos.y > ChainReaction::kShelfTop + ChainReaction::kDominoHalf.y);
    CHECK(!c.place(1, {8.05f, 4.0f}));  // on top of the first one
    CHECK(c.placed_count() == 1);
}

TEST(chain_reset_empties_the_gap_and_replays_identically) {
    ChainReaction a;
    for (int i = 0; i < 3; ++i) a.place(i, {7.25f + 0.75f * static_cast<Real>(i), 4.0f});
    const Status first = play(a);
    const double first_time = a.sim_time;
    a.reset();
    CHECK(a.placed_count() == 0);
    CHECK(a.phase == Phase::Setup);
    CHECK(a.world.bodies.size() == ChainReaction().world.bodies.size());
    for (int i = 0; i < 3; ++i) a.place(i, {7.25f + 0.75f * static_cast<Real>(i), 4.0f});
    CHECK(play(a) == first);
    CHECK(a.sim_time == first_time);
}

// ---- Level 5: Impossible Tower ----

TEST(tower_nothing_built_or_too_short_fails) {
    Tower empty;
    CHECK(play(empty) == Status::Failed);

    Tower low;
    CHECK(low.place(3, {10, 3}));
    CHECK(low.place(4, {10, 5}));
    CHECK_NEAR(low.height(), 2.4, 0.01);
    CHECK(play(low) == Status::Failed);
    CHECK(low.best_height == 0);
}

TEST(tower_a_sound_stack_wins_only_after_the_hold_time) {
    Tower t;
    CHECK(t.place(3, {10, 3}));
    CHECK(t.place(0, {10, 6}, kPi / 2));  // a plank on end
    CHECK(t.place(4, {10, 8}));
    CHECK(t.height() > Tower::kGoal);
    t.start();
    while (t.sim_time < Tower::kHoldTime - 0.05) {
        t.step();
        CHECK(t.status == Status::Playing);  // tall from the first frame, and still not accepted
    }
    for (int i = 0; i < 120 * 5 && t.status == Status::Playing; ++i) t.step();
    CHECK(t.status == Status::Won);
    CHECK(t.win_time >= Tower::kHoldTime);
    CHECK(t.height() > Tower::kGoal);
    CHECK_NEAR(t.best_height, t.height(), 1e-4);
}

// Tall when the run starts, on the floor two seconds later: height on one frame proves nothing.
TEST(tower_that_is_tall_but_topples_fails) {
    Tower t;
    const Vec2 spots[] = {{10, 3}, {10.9f, 5}, {11.7f, 6}, {12.6f, 7}, {13.4f, 8}, {13.9f, 9}, {14.2f, 10.4f}};
    const int order[] = {3, 4, 0, 1, 2, 6, 5};
    for (int i = 0; i < 7; ++i) CHECK(t.place(order[i], spots[i]));
    CHECK(t.height() > Tower::kGoal);
    t.start();
    t.step();
    CHECK(t.height() > Tower::kGoal);
    CHECK(t.status == Status::Playing);
    for (int i = 0; i < 120 * 30 && t.status == Status::Playing; ++i) t.step();
    CHECK(t.status == Status::Failed);
    CHECK(t.height() < Tower::kGoal);
    CHECK(t.held() < Tower::kHoldTime);
    CHECK(t.best_height == 0);
}

TEST(tower_best_height_survives_reset_and_only_grows) {
    Tower t;
    t.place(0, {10, 4}, kPi / 2);
    t.place(1, {10, 8}, kPi / 2);  // two planks on end: 4.8 m
    CHECK(play(t) == Status::Won);
    const Real best = t.best_height;
    CHECK(best > 4.7f);
    t.reset();
    CHECK(t.placed_count() == 0);
    CHECK(t.best_height == best);
    t.place(3, {10, 3});
    CHECK(play(t) == Status::Failed);
    CHECK(t.best_height == best);
}

TEST(tower_a_shorter_win_does_not_lower_the_best) {
    Tower t;
    const Vec2 spots[] = {{10, 3}, {10, 5}, {10, 6}, {10, 7}, {10, 8}, {10, 9}, {10, 10.5f}};
    const int order[] = {3, 4, 0, 1, 2, 6, 5};
    for (int i = 0; i < 7; ++i) CHECK(t.place(order[i], spots[i]));  // everything laid flat: about 5.4 m
    CHECK(play(t) == Status::Won);
    const Real best = t.best_height;
    CHECK(best > 5.3f);
    t.reset();
    t.place(0, {10, 4}, kPi / 2);
    t.place(1, {10, 8}, kPi / 2);  // 4.8 m: passes, but lower
    CHECK(play(t) == Status::Won);
    CHECK(t.height() < 5.0f);
    CHECK(t.best_height == best);
}

// ---- Level 2: Bridge Builder ----

namespace {

struct Edge {
    int c1, r1, c2, r2;
};
bool add(Bridge& b, std::initializer_list<Edge> edges) {
    bool ok = true;
    for (const Edge& e : edges) ok = b.add_beam(Bridge::node(e.c1, e.r1), Bridge::node(e.c2, e.r2)) && ok;
    return ok;
}
constexpr Edge kDeck[] = {{0, 1, 1, 1}, {1, 1, 2, 1}, {2, 1, 3, 1}, {3, 1, 4, 1}};
constexpr Edge kEnds[] = {{0, 0, 1, 1}, {4, 0, 4, 1}, {4, 0, 3, 1}};
void add_deck_and_ends(Bridge& b) {
    for (const Edge& e : kDeck) CHECK(b.add_beam(Bridge::node(e.c1, e.r1), Bridge::node(e.c2, e.r2)));
    for (const Edge& e : kEnds) CHECK(b.add_beam(Bridge::node(e.c1, e.r1), Bridge::node(e.c2, e.r2)));
}

}  // namespace

TEST(bridge_with_no_bridge_the_car_falls) {
    Bridge b;
    CHECK(play(b) == Status::Failed);
    CHECK(b.chassis().pos.y < Bridge::kFallY);
    CHECK(b.chassis().pos.x > Bridge::kLeftEdge);  // it drove off the edge under its own power
}

TEST(bridge_a_bare_roadway_does_not_hold) {
    Bridge b;
    for (const Edge& e : kDeck) CHECK(b.add_beam(Bridge::node(e.c1, e.r1), Bridge::node(e.c2, e.r2)));
    CHECK(play(b) == Status::Failed);

    Bridge braced_ends;  // seven beams: one freedom left, and it sags under the car
    add_deck_and_ends(braced_ends);
    CHECK(play(braced_ends) == Status::Failed);
    CHECK(braced_ends.chassis().pos.x < Bridge::kFinishX);
}

TEST(bridge_triangulated_designs_carry_the_car) {
    Bridge left;  // a triangle under the left half
    add_deck_and_ends(left);
    CHECK(add(left, {{0, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 2, 1}}));
    CHECK(play(left) == Status::Won);
    CHECK(left.chassis().pos.x >= Bridge::kFinishX);
    CHECK(left.win_time > 2.0);
    for (int i = 0; i < 240; ++i) left.step();
    CHECK(left.chassis().pos.x < Bridge::kFinishX + 1.0f);  // it braked rather than driving on into the wall
    CHECK(left.status == Status::Won);

    Bridge centre;  // a different answer: a V under the middle
    add_deck_and_ends(centre);
    CHECK(add(centre, {{1, 1, 2, 0}, {2, 0, 3, 1}, {2, 0, 2, 1}}));
    CHECK(play(centre) == Status::Won);
}

TEST(bridge_one_beam_short_of_rigid_fails) {
    Bridge b;  // the left-triangle design without its post on the right
    for (const Edge& e : kDeck) CHECK(b.add_beam(Bridge::node(e.c1, e.r1), Bridge::node(e.c2, e.r2)));
    CHECK(add(b, {{0, 0, 1, 1}, {4, 0, 3, 1}, {0, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 2, 1}}));
    CHECK(play(b) == Status::Failed);
}

TEST(bridge_beams_are_bodies_joined_by_real_hinges) {
    Bridge b;
    const size_t bodies = b.world.bodies.size(), joints = b.world.joints.size();
    CHECK(b.add_beam(Bridge::node(0, 1), Bridge::node(1, 1)));  // cliff corner to a free point
    CHECK(b.world.bodies.size() == bodies + 1);
    CHECK(b.world.joints.size() == joints + 1);                 // pinned to the world at the cliff
    CHECK(b.world.joints.back().type == JointType::Revolute);
    CHECK(b.world.joints.back().a == -1);
    CHECK(b.add_beam(Bridge::node(1, 1), Bridge::node(2, 1)));
    CHECK(b.world.joints.size() == joints + 2);                 // and to each other where they meet
    CHECK(b.world.joints.back().a >= 0 && b.world.joints.back().b >= 0);
    const Body& beam = b.world.bodies.back();
    CHECK(beam.type == BodyType::Dynamic);
    CHECK_NEAR(beam.pos.x, 9.0, 1e-4);
    CHECK_NEAR(beam.pos.y, Bridge::kRoadY, 1e-4);
}

TEST(bridge_connection_rules) {
    CHECK(Bridge::is_anchor(Bridge::node(0, 0)) && Bridge::is_anchor(Bridge::node(0, 1)) && Bridge::is_anchor(Bridge::node(4, 0)));
    CHECK(!Bridge::is_anchor(Bridge::node(4, 1)));  // the crumbling edge
    CHECK(!Bridge::is_anchor(Bridge::node(0, 2)) && !Bridge::is_anchor(Bridge::node(2, 1)));
    CHECK(Bridge::can_connect(Bridge::node(1, 1), Bridge::node(2, 1)));
    CHECK(Bridge::can_connect(Bridge::node(1, 1), Bridge::node(2, 0)));   // one diagonal
    CHECK(!Bridge::can_connect(Bridge::node(1, 1), Bridge::node(3, 1)));  // two apart
    CHECK(!Bridge::can_connect(Bridge::node(1, 0), Bridge::node(1, 2)));
    CHECK(!Bridge::can_connect(Bridge::node(0, 0), Bridge::node(0, 1)));  // both ends in rock
    CHECK(!Bridge::can_connect(Bridge::node(2, 1), Bridge::node(2, 1)));
    CHECK(!Bridge::can_connect(Bridge::node(2, 1), -1));
    CHECK(!Bridge::can_connect(Bridge::node(2, 1), 99));

    Bridge b;
    CHECK(b.add_beam(Bridge::node(1, 1), Bridge::node(2, 1)));
    CHECK(!b.add_beam(Bridge::node(2, 1), Bridge::node(1, 1)));  // already there, whichever way round
    CHECK(!b.add_beam(Bridge::node(1, 1), Bridge::node(3, 1)));
    CHECK(b.beams.size() == 1);
    CHECK(b.remove_beam(Bridge::node(2, 1), Bridge::node(1, 1)));
    CHECK(!b.remove_beam(Bridge::node(2, 1), Bridge::node(1, 1)));
    CHECK(b.beams.empty());
    CHECK(b.world.bodies.size() == Bridge().world.bodies.size());
}

TEST(bridge_beam_limit) {
    Bridge b;
    int added = 0;
    for (int c = 0; c < Bridge::kColumns - 1; ++c)
        for (int r = 0; r < Bridge::kRows; ++r) {
            added += b.add_beam(Bridge::node(c, r), Bridge::node(c + 1, r)) ? 1 : 0;
            if (r + 1 < Bridge::kRows) added += b.add_beam(Bridge::node(c, r), Bridge::node(c + 1, r + 1)) ? 1 : 0;
        }
    CHECK(added == Bridge::kMaxBeams);
    CHECK(static_cast<int>(b.beams.size()) == Bridge::kMaxBeams);
    CHECK(!b.add_beam(Bridge::node(3, 2), Bridge::node(4, 2)));
    CHECK(b.remove_beam(b.beams[0].a, b.beams[0].b));
    CHECK(b.add_beam(Bridge::node(3, 2), Bridge::node(4, 2)));  // room again
}

TEST(bridge_is_built_with_the_pointer_and_frozen_during_a_run) {
    Bridge b;
    const Vec2 p = Bridge::node_position(Bridge::node(1, 1)), q = Bridge::node_position(Bridge::node(2, 1));
    b.press(p + Vec2{0.2f, -0.1f});  // near enough
    b.release(q + Vec2{-0.1f, 0.2f});
    CHECK(b.beams.size() == 1);
    b.press(p);
    b.release({9.0f, 8.5f});  // let go over nothing
    CHECK(b.beams.size() == 1);
    b.press(p);
    CHECK(b.cancel());
    b.release(Bridge::node_position(Bridge::node(1, 2)));
    CHECK(b.beams.size() == 1);

    b.start();
    b.press((p + q) * 0.5f);
    b.release((p + q) * 0.5f);
    CHECK(b.beams.size() == 1);
    CHECK(!b.add_beam(Bridge::node(2, 1), Bridge::node(3, 1)));
    CHECK(!b.remove_beam(Bridge::node(1, 1), Bridge::node(2, 1)));
    b.stop();

    b.press((p + q) * 0.5f);  // a click on a beam takes it away
    b.release((p + q) * 0.5f);
    CHECK(b.beams.empty());
}

TEST(bridge_reset_and_replay) {
    Bridge a;
    add_deck_and_ends(a);
    add(a, {{0, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 2, 1}});
    CHECK(play(a) == Status::Won);
    const double t = a.win_time;
    a.stop();
    CHECK(a.beams.size() == 10);  // the design is kept
    CHECK(a.chassis().pos.x < Bridge::kLeftEdge);  // and the car is back at the start
    CHECK(play(a) == Status::Won);
    CHECK(a.win_time == t);
    a.reset();
    CHECK(a.beams.empty());
    CHECK(a.world.bodies.size() == Bridge().world.bodies.size());
    CHECK(a.world.joints.size() == Bridge().world.joints.size());
}
