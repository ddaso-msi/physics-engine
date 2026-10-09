#pragma once
// Level 1, Knock It Down. A pyramid of blocks stands on a raised platform; the player has a few balls
// to knock every block off it. A block counts as down once its centre is below the platform's top,
// which it cannot be while it is still resting on the platform.
#include "level.hpp"

#include <vector>

namespace puzzle {

struct KnockDown {
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

    World world;
    std::vector<int> targets;      // the blocks, as indices into world.bodies
    size_t first_ball = 0;         // bodies from here on are balls the player fired
    int shots_left = kShots;
    int fallen = 0;                // targets below the fall line after the latest step
    double sim_time = 0;           // simulated seconds since reset
    double last_shot_time = 0;
    double win_time = 0;           // sim_time when the last block went down
    Status status = Status::Playing;

    KnockDown() { reset(); }

    // Back to the initial state: a new world, full shots, zero time.
    void reset();
    // One fixed step of physics, then the objective.
    void step(Real dt = kTimeStep);

    // The velocity a ball would get if the player aimed at world point `aim`.
    Vec2 launch_velocity(Vec2 aim) const;
    bool can_fire() const;
    // Fires at `aim` if allowed; returns whether a ball was launched.
    bool fire(Vec2 aim);
    int shots_used() const { return kShots - shots_left; }
};

}  // namespace puzzle
