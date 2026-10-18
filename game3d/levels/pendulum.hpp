#pragma once
// Level 3, Pendulum Smash. A heavy weight hangs on a rod from a fixed pivot; a ball sits on a table in its
// path and a bin waits beyond the table. The player pulls the ball back to any point on the sphere the rod
// sweeps: how far back sets the speed, and which way sets the direction. Too gentle and the ball drops off
// the table short of the bin, too hard and it flies over, a few degrees off line and it lands beside it.
// (The target is a ball, not a crate, because a ball goes off along the line between the two centres; a crate
// hit anywhere on its face goes straight ahead whatever the aim.)
#include "level.hpp"
#include "objectives.hpp"

namespace puzzle3d {

class Pendulum : public Level {
public:
    static constexpr Vec3 kPivot{0, 7.5f, 0};
    static constexpr Real kRod = 5.0f, kBobRadius = 0.5f;
    // The pull: `tilt` is the rod's angle from straight down (radians), `heading` the direction it points,
    // measured about the vertical from -x (straight back from the bin) toward +z.
    static constexpr Real kMinTilt = 0.08f, kMaxTilt = 1.57f, kStartTilt = 0.35f;
    static constexpr Vec3 kZoneLo{7.0f, kGroundY, -1.1f}, kZoneHi{8.8f, 1.5f, 1.1f};
    static constexpr double kRestTime = 0.75;
    static constexpr Real kTargetRadius = 0.35f;
    // Where the weight's centre is level with the ball's, 25.5 degrees up the far side, the two just touch.
    static constexpr Vec3 kTargetStart{3.0f, 2.987f, 0};
    static constexpr Real kTableTop = kTargetStart.y - kTargetRadius - 0.001f;

    Real tilt = kStartTilt, heading = 0;
    int bob = -1, target = -1;

    Pendulum() { reset(); }

    const char* name() const override { return "Pendulum Smash"; }
    const char* objective() const override { return "Swing the weight to knock the ball into the bin."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    View view() const override { return {{3.0f, 3.6f, 0}, -0.32f, 0.3f, 16.5f}; }
    Role role(size_t body) const override;
    void overlay(Overlay& out, const Ray* pointer) const override;
    void press(const Ray& ray) override;
    void drag(const Ray& ray) override;
    void release(const Ray&) override { dragging_ = false; }
    void rotate(int direction) override;  // the heading, one degree at a time
    bool cancel() override;

    static Vec3 bob_position(Real tilt, Real heading);
    // Sets the pull (clamped to the allowed tilts). Only during setup.
    void set_pull(Real new_tilt, Real new_heading);
    bool target_in_bin() const { return inside(world, target, kZoneLo, kZoneHi); }

protected:
    void clear_setup() override {
        tilt = kStartTilt;
        heading = 0;
    }
    void build() override;
    void evaluate() override;

private:
    bool dragging_ = false;
    Real tilt_before_ = kStartTilt, heading_before_ = 0;
    Hold resting_, stalled_;
};

}  // namespace puzzle3d
