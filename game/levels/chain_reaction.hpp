#pragma once
// Level 4, Chain Reaction. A ball rolls down a ramp into a row of dominoes; the last domino nudges a
// second ball off the end of the shelf into a bin. Some dominoes are missing from the middle of the row
// and the player has a few spares to bridge the gap. Nothing is scripted: if the spares are too far
// apart, or too close, the chain simply stops.
#include "objectives.hpp"
#include "pieces.hpp"

namespace puzzle {

class ChainReaction : public PieceLevel {
public:
    static constexpr Real kShelfTop = 3.0f, kShelfLeft = 3.6f, kShelfRight = 12.8f;
    static constexpr Vec2 kDominoHalf{0.075f, 0.5f};
    static constexpr Real kGapLeft = 6.5f, kGapRight = 9.5f;  // the dominoes either side of the gap stand here
    static constexpr Vec2 kZoneLo{13.0f, kFloorTop}, kZoneHi{15.8f, 2.0f};
    static constexpr double kRestTime = 1.0;  // the ball must stay in the bin this long

    int starter = -1, target = -1;  // the two balls, as indices into world.bodies

    ChainReaction() { reset(); }

    const char* name() const override { return "Chain Reaction"; }
    const char* objective() const override { return "Bridge the gap so the last ball lands in the bin."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    void overlay(Overlay& out, Vec2 pointer, bool pointer_active) const override;

    bool target_in_bin() const { return inside(world, target, kZoneLo, kZoneHi); }

protected:
    std::vector<Piece> initial_pieces() const override;
    void build_scenery() override;
    phys::AABB build_area() const override { return {{kGapLeft + 0.3f, kShelfTop}, {kGapRight - 0.3f, kShelfTop + 2.2f}}; }
    Role scenery_role(size_t body) const override;
    void evaluate() override;

private:
    Hold resting_, stalled_;
};

}  // namespace puzzle
