#include "levels/knock_down.hpp"

#include "objectives.hpp"

#include <scenes.hpp>

#include <algorithm>
#include <cstdio>

namespace puzzle {

std::span<const char* const> KnockDown::help() const {
    static constexpr const char* kLines[] = {
        "Move the mouse to aim. The dots show the path.",
        "Aim further from the launcher for more speed.",
        "Click to fire. You have three balls.",
        "Every block must drop below the fall line.",
    };
    return kLines;
}

const char* KnockDown::hint() const {
    if (status == Status::Won) return "Solved.";
    if (status == Status::Failed) return "Out of balls. R to try again.";
    if (shots_left == 0) return "Last ball away. Watch what it does.";
    return can_fire() ? "Aim with the mouse, click to fire." : "Reloading...";
}

std::string KnockDown::stats() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "Balls %d/%d   Down %d/%d", shots_left, kShots, fallen, static_cast<int>(targets.size()));
    return buf;
}

std::string KnockDown::result() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%d ball%s, %.1f s", shots_used(), shots_used() == 1 ? "" : "s", win_time);
    return buf;
}

Role KnockDown::role(size_t body) const {
    if (body >= first_ball) return Role::Tool;
    if (std::find(targets.begin(), targets.end(), static_cast<int>(body)) == targets.end()) return Role::Scenery;
    return world.bodies[body].pos.y < kFallLine ? Role::Done : Role::Target;
}

void KnockDown::overlay(Overlay& out, Vec2 pointer, bool pointer_active) const {
    out.lines.push_back({{kPlatformCenter.x - 2.4f, kFallLine}, {kPlatformCenter.x + 2.4f, kFallLine}, Role::Target, true});
    out.labels.push_back({{kPlatformCenter.x + 2.5f, kFallLine + 0.12f}, "fall line", Role::Target});

    // The launcher, with the next ball loaded while there is one.
    out.dots.push_back({kMuzzle, 0.36f, Role::Scenery});
    if (shots_left > 0 && status == Status::Playing) out.dots.push_back({kMuzzle, kBallRadius, can_fire() ? Role::Tool : Role::Guide});
    if (!pointer_active || !can_fire()) return;

    // The same arithmetic the engine uses for a body in free flight, so the dots are where the ball goes.
    const Vec2 v0 = launch_velocity(pointer);
    Vec2 p = kMuzzle, v = v0;
    for (int i = 1; i <= 150; ++i) {
        v += world.gravity * kTimeStep;
        p += v * kTimeStep;
        if (p.y < kFloorTop || p.x < 0 || p.x > kWorldW) break;
        if (i % 5 == 0) out.dots.push_back({p, 0.04f, Role::Guide});
    }
    out.lines.push_back({kMuzzle, kMuzzle + v0.normalized() * 0.75f, Role::Tool, false});
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.0f m/s", static_cast<double>(v0.length()));
    out.labels.push_back({kMuzzle + Vec2{-0.7f, 0.9f}, buf, Role::Tool});
}

void KnockDown::clear_setup() {
    shots_left = kShots;
    last_shot_time = 0;
}

void KnockDown::build() {
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
    fallen = 0;
}

void KnockDown::evaluate() {
    fallen = count_below(world, targets, kFallLine);
    if (fallen == static_cast<int>(targets.size())) {  // a late block can still turn a failure into a win
        win();
    } else if (shots_left == 0) {
        const double since = sim_time - last_shot_time;
        const bool at_rest = world.stats().awake_bodies == 0 && since >= 0.5;
        if (at_rest || since >= static_cast<double>(kSettleTime)) fail();
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
    mark_tried();
    return true;
}

}  // namespace puzzle
