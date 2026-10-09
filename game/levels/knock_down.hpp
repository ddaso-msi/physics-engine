#pragma once
// Level 1, Knock It Down. A pyramid of blocks stands on a raised platform; the player has a few balls
// to knock every block off it. A block counts as down once its centre is below the platform's top,
// which it cannot be while it is still resting on the platform.
//
// There is no setup phase: physics runs the whole time and the player fires into it.
#include "level.hpp"

#include <vector>

namespace puzzle {

class KnockDown : public Level {
public:
    static constexpr int kShots = 3;
    static constexpr Vec2 kMuzzle{2.0f, 2.6f};  // where a ball starts
    static constexpr Real kBallRadius = 0.25f;
    // Launch speed grows with how far the aim point is from the muzzle, between these limits. The top
    // speed moves the ball 17 cm a step, well under a block's width: the engine's continuous collision
    // only covers static bodies, so a faster ball could pass through a block.
    static constexpr Real kMinSpeed = 4.0f, kMaxSpeed = 20.0f, kSpeedPerMeter = 2.2f;
    static constexpr Real kReloadTime = 0.4f;   // simulated seconds between shots: the last ball is clear by then
    static constexpr Real kSettleTime = 6.0f;   // after the last shot, how long the blocks get to fall
    static constexpr Vec2 kPlatformCenter{14.5f, 3.3f}, kPlatformHalf{1.3f, 0.2f};
    static constexpr Real kFallLine = kPlatformCenter.y + kPlatformHalf.y - 0.1f;

    std::vector<int> targets;      // the blocks, as indices into world.bodies
    size_t first_ball = 0;         // bodies from here on are balls the player fired
    int shots_left = kShots;
    int fallen = 0;                // targets below the fall line after the latest step
    double last_shot_time = 0;

    KnockDown() { reset(); }

    const char* name() const override { return "Knock It Down"; }
    const char* objective() const override { return "Knock every block off the platform."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    bool has_setup() const override { return false; }
    Role role(size_t body) const override;
    void overlay(Overlay& out, Vec2 pointer, bool pointer_active) const override;
    void press(Vec2 p) override { fire(p); }

    // The velocity a ball would get if the player aimed at world point `aim`.
    Vec2 launch_velocity(Vec2 aim) const;
    bool can_fire() const;
    // Fires at `aim` if allowed; returns whether a ball was launched.
    bool fire(Vec2 aim);
    int shots_used() const { return kShots - shots_left; }

protected:
    void clear_setup() override;
    void build() override;
    void evaluate() override;
};

}  // namespace puzzle
