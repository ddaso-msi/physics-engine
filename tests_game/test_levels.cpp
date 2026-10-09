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
