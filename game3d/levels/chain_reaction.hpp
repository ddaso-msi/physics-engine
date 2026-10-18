#pragma once
// Level 4, Chain Reaction. Dominoes stand along a half circle. A ball rolls down a ramp into the first;
// the last tips a second ball into a bin. Three dominoes are missing from the middle of the curve and the
// player has three spares: they have to be spaced along the curve AND turned to follow it, because a domino
// falls the way it faces.
#include "objectives.hpp"
#include "pieces.hpp"

namespace puzzle3d {

class ChainReaction : public PieceLevel {
public:
    static constexpr Real kRadius = 3.0f;     // the curve: a half circle about the origin, from +x round through +z to -x
    static constexpr Real kPitch = 0.65f;     // metres of curve between neighbouring dominoes
    static constexpr int kSlots = 14;         // positions along it...
    static constexpr int kGapFirst = 2, kGapCount = 3;  // ...of which these are left empty
    static constexpr Vec3 kDominoHalf{0.075f, 0.5f, 0.3f};  // thin along the curve, tall, wide across it
    static constexpr Vec3 kZoneLo{-3.9f, kGroundY, -2.9f}, kZoneHi{-2.1f, 1.5f, -1.1f};
    static constexpr double kRestTime = 1.0;

    int starter = -1, target = -1;

    ChainReaction() { reset(); }

    const char* name() const override { return "Chain Reaction"; }
    const char* objective() const override { return "Fill the gap so the last ball lands in the bin."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    View view() const override { return {{0, 0.5f, 0.6f}, 0.35f, 0.85f, 12.0f}; }
    void overlay(Overlay& out, const Ray* pointer) const override;

    // Slot k on the curve: where a domino there stands, and the turn about the vertical that faces it along
    // the curve.
    static Vec3 slot_position(int k);
    static Real slot_yaw(int k);
    bool target_in_bin() const { return inside(world, target, kZoneLo, kZoneHi); }

protected:
    std::vector<Piece> initial_pieces() const override;
    void build_scenery() override;
    phys3d::AABB build_area() const override;
    Role scenery_role(size_t body) const override;
    void evaluate() override;

private:
    std::vector<int> chain_;  // the dominoes and the last ball: what has to keep moving for the chain to be alive
    Hold resting_, stalled_;
};

}  // namespace puzzle3d
