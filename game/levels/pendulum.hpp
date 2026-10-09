#pragma once
// Level 3, Pendulum Smash. A heavy weight hangs on a rod from a fixed pivot; a crate sits on a table in
// its path. The player chooses how far back to pull the weight. Too little and the crate stops on the
// table or drops short; too much and it flies past the bin.
#include "level.hpp"
#include "objectives.hpp"

namespace puzzle {

class Pendulum : public Level {
public:
    static constexpr Vec2 kPivot{6.6f, 9.6f};
    static constexpr Real kRod = 5.2f, kBobRadius = 0.5f;
    // Release angle from straight down, in radians; negative is to the left, the only side offered.
    static constexpr Real kMinAngle = -1.57f, kMaxAngle = -0.1f, kStartAngle = -0.35f;
    static constexpr Vec2 kZoneLo{12.7f, kFloorTop}, kZoneHi{15.1f, 2.0f};
    static constexpr double kRestTime = 0.75;  // the crate must stay in the bin, nearly still, this long

    Real angle = kStartAngle;
    int bob = -1, crate = -1;  // indices into world.bodies

    Pendulum() { reset(); }

    const char* name() const override { return "Pendulum Smash"; }
    const char* objective() const override { return "Swing the weight to knock the crate into the bin."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    Role role(size_t body) const override;
    void overlay(Overlay& out, Vec2 pointer, bool pointer_active) const override;
    void press(Vec2 p) override;
    void drag(Vec2 p) override;
    void release(Vec2) override { dragging_ = false; }
    void rotate(int direction) override;  // one degree at a time, for fine tuning
    bool cancel() override;

    // Sets the release angle (clamped to the allowed range). Only during setup.
    void set_angle(Real a);
    bool crate_in_bin() const { return inside(world, crate, kZoneLo, kZoneHi); }

protected:
    void clear_setup() override { angle = kStartAngle; }
    void build() override;
    void evaluate() override;

private:
    bool dragging_ = false;
    Real angle_before_drag_ = kStartAngle;
    Hold resting_, stalled_;
};

}  // namespace puzzle
