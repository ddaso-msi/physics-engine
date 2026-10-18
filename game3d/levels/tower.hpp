#pragma once
// Level 5, Impossible Tower. Build the tallest tower you can from a small set of mismatched pieces on a
// square pad. It only counts if it is still above the line after standing for a while with physics on:
// height is measured from the simulation, continuously, never from a single frame.
#include "objectives.hpp"
#include "pieces.hpp"

namespace puzzle3d {

class Tower : public PieceLevel {
public:
    static constexpr Real kGoal = 4.5f;       // metres above the ground
    static constexpr double kHoldTime = 5.0;  // simulated seconds the tower must stay above the line
    static constexpr Real kPad = 2.0f;        // half the width of the square it has to be built on

    Real best_height = 0;  // the tallest tower that passed, across attempts (reset() keeps it)

    Tower() { reset(); }

    const char* name() const override { return "Impossible Tower"; }
    const char* objective() const override { return "Build a tower that stays above the line."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    View view() const override { return {{-1.5f, 2.4f, 0}, 0.55f, 0.32f, 15.0f}; }
    void overlay(Overlay& out, const Ray* pointer) const override;

    Real height() const;  // of the tallest point of any placed piece above the ground, now
    double held() const { return standing_.held; }

protected:
    std::vector<Piece> initial_pieces() const override;
    void build_scenery() override;
    phys3d::AABB build_area() const override { return {{-kPad, kGroundY, -kPad}, {kPad, 12.0f, kPad}}; }
    void evaluate() override;

private:
    Hold standing_, collapsed_;
};

}  // namespace puzzle3d
