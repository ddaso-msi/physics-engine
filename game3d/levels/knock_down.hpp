#pragma once
// Level 1, Knock It Down. A pyramid of crates (three layers: 9, 4 and 1) stands on a raised platform; the
// player has a few balls to knock every crate off it. A crate counts as down once its centre is below the
// platform's top, which it cannot be while it is still resting on the platform.
//
// There is no setup phase: physics runs the whole time and the player fires into it.
#include "level.hpp"

#include <vector>

namespace puzzle3d {

class KnockDown : public Level {
public:
    static constexpr int kShots = 4;
    static constexpr Vec3 kMuzzle{0, 2.0f, 0};
    static constexpr Real kBallRadius = 0.3f;
    // Launch speed grows with how far the aim point is from the muzzle, between these limits. The top speed
    // moves a ball 18 cm a step, well under a crate's width: continuous collision only covers static bodies,
    // so a faster ball could pass through a crate.
    static constexpr Real kMinSpeed = 4.0f, kMaxSpeed = 22.0f, kSpeedPerMeter = 1.8f;
    static constexpr Real kReloadTime = 0.4f;
    static constexpr Real kSettleTime = 6.0f;  // after the last shot, how long the crates get to fall
    static constexpr Vec3 kPlatformCenter{12.0f, 2.35f, 0}, kPlatformHalf{1.4f, 0.15f, 1.4f};
    static constexpr Real kFallLine = kPlatformCenter.y + kPlatformHalf.y - 0.1f;

    std::vector<int> targets;
    size_t first_ball = 0;
    int shots_left = kShots;
    int fallen = 0;
    double last_shot_time = 0;

    KnockDown() { reset(); }

    const char* name() const override { return "Knock It Down"; }
    const char* objective() const override { return "Knock every crate off the platform."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    View view() const override { return {{7.5f, 2.4f, 0}, -1.05f, 0.33f, 12.5f}; }
    bool has_setup() const override { return false; }
    Role role(size_t body) const override;
    void overlay(Overlay& out, const Ray* pointer) const override;
    void press(const Ray& ray) override { fire(aim_point(ray)); }

    // Where a ball would be aimed for the pointer on this ray: the first thing it meets, or 20 m out.
    Vec3 aim_point(const Ray& ray) const;
    Vec3 launch_velocity(Vec3 aim) const;
    bool can_fire() const;
    bool fire(Vec3 aim);
    int shots_used() const { return kShots - shots_left; }

protected:
    void clear_setup() override;
    void build() override;
    void evaluate() override;
};

}  // namespace puzzle3d
