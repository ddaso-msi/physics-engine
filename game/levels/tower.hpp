#pragma once
// Level 5, Impossible Tower. Build the tallest tower you can from a small set of mismatched pieces. It
// only counts if it is still above the line after standing for a while with physics switched on:
// height is measured from the simulation, continuously, never from a single frame.
#include "objectives.hpp"
#include "pieces.hpp"

namespace puzzle {

class Tower : public PieceLevel {
public:
    static constexpr Real kGoal = 4.5f;        // metres above the floor
    static constexpr double kHoldTime = 5.0;   // simulated seconds the tower must stay above the line

    Real best_height = 0;  // the tallest tower that passed, across attempts (reset() keeps it)

    Tower() { reset(); }

    const char* name() const override { return "Impossible Tower"; }
    const char* objective() const override { return "Build a tower that stays above the line."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    void overlay(Overlay& out, Vec2 pointer, bool pointer_active) const override;

    // Height of the tallest placed piece above the floor, right now.
    Real height() const;
    // How long the tower has been above the line without a break, in the current run.
    double held() const { return standing_.held; }

protected:
    std::vector<Piece> initial_pieces() const override;
    void build_scenery() override;
    phys::AABB build_area() const override { return {{5.6f, kFloorTop}, {kWorldW - 0.6f, kWorldH - 0.1f}}; }
    Real rotation_step() const override { return phys::kPi / 2; }
    void evaluate() override;

private:
    Hold standing_, collapsed_;
};

}  // namespace puzzle
