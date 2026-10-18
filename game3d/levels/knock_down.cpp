#include "levels/knock_down.hpp"

#include "objectives.hpp"

#include <algorithm>
#include <cstdio>

namespace puzzle3d {

std::span<const char* const> KnockDown::help() const {
    static constexpr const char* kLines[] = {
        "Point the mouse to aim. The dots show the path.",
        "Aim further from the launcher for more speed.",
        "Click to fire. You have four balls.",
        "Every crate must drop below the platform's top.",
        "Right-drag turns the camera; the wheel zooms.",
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
    if (std::find(targets.begin(), targets.end(), static_cast<int>(body)) == targets.end()) return Level::role(body);
    return world.bodies[body].pos.y < kFallLine ? Role::Done : Role::Target;
}

void KnockDown::overlay(Overlay& out, const Ray* pointer) const {
    // The fall line: a dashed square just under the platform's top.
    const Real h = kPlatformHalf.x + 0.25f, y = kFallLine;
    const Vec3 c{kPlatformCenter.x, y, kPlatformCenter.z};
    const Vec3 corners[4] = {c + Vec3{-h, 0, -h}, c + Vec3{h, 0, -h}, c + Vec3{h, 0, h}, c + Vec3{-h, 0, h}};
    for (int i = 0; i < 4; ++i) out.lines.push_back({corners[i], corners[(i + 1) % 4], Role::Target, true});
    out.labels.push_back({c + Vec3{0, -0.3f, h + 0.2f}, "fall line", Role::Target});

    if (shots_left > 0 && status == Status::Playing) {
        Body loaded = Body::solid_sphere(kBallRadius, 1, kMuzzle);
        out.ghosts.push_back({loaded, can_fire() ? Role::Tool : Role::Guide});
    }
    if (!pointer || !can_fire()) return;
    // The same arithmetic the engine uses for a body in free flight, so the dots are where the ball goes.
    const Vec3 v0 = launch_velocity(aim_point(*pointer));
    Vec3 p = kMuzzle, v = v0;
    for (int i = 1; i <= 180; ++i) {
        v += world.gravity * kTimeStep;
        p += v * kTimeStep;
        if (p.y < kGroundY) break;
        if (i % 5 == 0) out.dots.push_back({p, 0.05f, Role::Guide});
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.0f m/s", static_cast<double>(v0.length()));
    out.labels.push_back({kMuzzle + Vec3{0, 0.9f, 0}, buf, Role::Tool});
}

Vec3 KnockDown::aim_point(const Ray& ray) const {
    if (const auto hit = raycast(world, ray)) return hit->point;
    return ray.at(20.0f);
}

void KnockDown::clear_setup() {
    shots_left = kShots;
    last_shot_time = 0;
}

void KnockDown::build() {
    targets.clear();
    add_ground(world);
    world.add(Body::fixed_box({0.45f, 0.8f, 0.45f}, {kMuzzle.x, 0.8f, kMuzzle.z}));  // the launcher's pedestal
    world.add(Body::fixed_box(kPlatformHalf, kPlatformCenter));
    const Real pillar = (kPlatformCenter.y - kPlatformHalf.y - kGroundY) * 0.5f;
    world.add(Body::fixed_box({0.3f, pillar, 0.3f}, {kPlatformCenter.x, kGroundY + pillar, kPlatformCenter.z}));

    // Three layers, 3 x 3, 2 x 2 and 1, a millimetre apart so nothing starts out overlapping.
    constexpr Real half = 0.35f, pitch = 0.75f;
    const Real top = kPlatformCenter.y + kPlatformHalf.y;
    for (int layer = 0; layer < 3; ++layer) {
        const int n = 3 - layer;
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k) {
                const Real off = static_cast<Real>(n - 1) * 0.5f;
                const Vec3 p{kPlatformCenter.x + (static_cast<Real>(i) - off) * pitch, top + half + static_cast<Real>(layer) * (2 * half + 0.001f),
                             kPlatformCenter.z + (static_cast<Real>(k) - off) * pitch};
                Body crate = Body::solid_box({half, half, half}, 1, p);
                crate.restitution = 0.1f;
                crate.friction = 0.5f;
                targets.push_back(world.add(crate));
            }
    }
    first_ball = world.bodies.size();
    fallen = 0;
}

void KnockDown::evaluate() {
    fallen = count_below(world, targets, kFallLine);
    if (fallen == static_cast<int>(targets.size())) {
        win();
    } else if (shots_left == 0) {
        const double since = sim_time - last_shot_time;
        const bool at_rest = world.stats().awake_bodies == 0 && since >= 0.5;
        if (at_rest || since >= static_cast<double>(kSettleTime)) fail();
    }
}

Vec3 KnockDown::launch_velocity(Vec3 aim) const {
    const Vec3 d = aim - kMuzzle;
    const Real len = d.length();
    if (len < 1e-4f) return {kMinSpeed, 0, 0};
    return d * (std::clamp(len * kSpeedPerMeter, kMinSpeed, kMaxSpeed) / len);
}

bool KnockDown::can_fire() const {
    if (status != Status::Playing || shots_left <= 0) return false;
    return shots_left == kShots || sim_time - last_shot_time >= static_cast<double>(kReloadTime);
}

bool KnockDown::fire(Vec3 aim) {
    if (!can_fire()) return false;
    Body ball = Body::solid_sphere(kBallRadius, 6.0f, kMuzzle);
    ball.restitution = 0.2f;
    ball.friction = 0.5f;
    ball.vel = launch_velocity(aim);
    world.add(ball);
    --shots_left;
    last_shot_time = sim_time;
    mark_tried();
    return true;
}

}  // namespace puzzle3d
