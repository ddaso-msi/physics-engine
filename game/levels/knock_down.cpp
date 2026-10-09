#include "levels/knock_down.hpp"

#include "objectives.hpp"

#include <scenes.hpp>

#include <algorithm>

namespace puzzle {

void KnockDown::reset() {
    // A new World rather than truncate(), so nothing at all carries over from the previous attempt
    // (broad-phase state, remembered contacts, sleep timers).
    world = World{};
    targets.clear();
    add_room(world);
    using phys::scenes::add_static_box;
    add_static_box(world, {0.45f, 0.8f}, {kMuzzle.x, kFloorTop + 0.8f});  // the launcher's pedestal
    add_static_box(world, kPlatformHalf, kPlatformCenter);
    const Real pillar = (kPlatformCenter.y - kPlatformHalf.y - kFloorTop) * 0.5f;
    add_static_box(world, {0.2f, pillar}, {kPlatformCenter.x, kFloorTop + pillar});

    // A 3-2-1 pyramid. Rows are a millimetre apart so nothing starts out overlapping.
    constexpr Real half = 0.35f, pitch = 0.75f;
    const Real top = kPlatformCenter.y + kPlatformHalf.y;
    for (int row = 0; row < 3; ++row)
        for (int i = 0; i < 3 - row; ++i) {
            const Real x = kPlatformCenter.x + (static_cast<Real>(i) - static_cast<Real>(2 - row) * 0.5f) * pitch;
            const Real y = top + half + static_cast<Real>(row) * (2 * half + 0.001f);
            Body block(phys::Shape::make_polygon(phys::Polygon::box(half, half)), {x, y}, 0);
            block.restitution = 0.1f;
            block.friction = 0.5f;
            targets.push_back(world.add(block));
        }

    first_ball = world.bodies.size();
    shots_left = kShots;
    fallen = 0;
    sim_time = last_shot_time = win_time = 0;
    status = Status::Playing;
}

void KnockDown::step(Real dt) {
    world.step(dt);
    sim_time += static_cast<double>(dt);
    fallen = count_below(world, targets, kFallLine);
    if (status == Status::Won) return;
    if (fallen == static_cast<int>(targets.size())) {  // a late block can still turn a failure into a win
        status = Status::Won;
        win_time = sim_time;
    } else if (status == Status::Playing && shots_left == 0) {
        const double since = sim_time - last_shot_time;
        const bool at_rest = world.stats().awake_bodies == 0 && since >= 0.5;
        if (at_rest || since >= static_cast<double>(kSettleTime)) status = Status::Failed;
    }
}

Vec2 KnockDown::launch_velocity(Vec2 aim) const {
    const Vec2 d = aim - kMuzzle;
    const Real len = d.length();
    if (len < 1e-4f) return {kMinSpeed, 0};
    return d * (std::clamp(len * kSpeedPerMeter, kMinSpeed, kMaxSpeed) / len);
}

bool KnockDown::can_fire() const {
    if (status != Status::Playing || shots_left <= 0) return false;
    return shots_left == kShots || sim_time - last_shot_time >= static_cast<double>(kReloadTime);
}

bool KnockDown::fire(Vec2 aim) {
    if (!can_fire()) return false;
    Body ball(phys::Shape::make_circle(kBallRadius), kMuzzle, 0, phys::BodyType::Dynamic, 6.0f);
    ball.restitution = 0.2f;
    ball.friction = 0.5f;
    ball.vel = launch_velocity(aim);
    world.add(ball);
    --shots_left;
    last_shot_time = sim_time;
    return true;
}

}  // namespace puzzle
