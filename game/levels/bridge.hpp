#pragma once
// Level 2, Bridge Builder. A car has to cross a gap between two cliffs. The player joins a grid of
// connection points with a limited number of beams; every beam end is a real hinge, to the cliff where
// the point is on one, otherwise to the other beams that meet there. Beams along the roadway carry
// the car; the rest are bracing and the car passes in front of them.
//
// The engine has hinges but no welds, so a beam is only held by what it is pinned to, and only
// triangles keep their shape. The right cliff offers no fixed point at road level, so a plain row of
// beams hangs loose there: the roadway has to be braced back to the three points that do hold. Counting
// freedoms, the smallest bridge that stands is ten beams; the limit leaves two to spare.
#include "level.hpp"

#include <scenes.hpp>

#include <vector>

namespace puzzle {

class Bridge : public Level {
public:
    static constexpr int kColumns = 5, kRows = 3;  // connection points; row 1 is the roadway
    static constexpr Real kLeftEdge = 6.0f, kSpacing = 2.0f, kRoadY = 3.9f;
    static constexpr Real kCliffTop = 4.0f;
    static constexpr int kMaxBeams = 12;
    static constexpr Real kFinishX = 16.5f;  // the car's body has to get this far along the right cliff
    static constexpr Real kFallY = 2.6f;     // ...without ever dropping this low

    struct Beam {
        int a, b;  // connection point indices, a < b
        bool operator==(const Beam&) const = default;
    };
    std::vector<Beam> beams;
    phys::scenes::Car car;

    Bridge() { reset(); }

    const char* name() const override { return "Bridge Builder"; }
    const char* objective() const override { return "Get the car across the gap."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    Role role(size_t body) const override;
    void overlay(Overlay& out, Vec2 pointer, bool pointer_active) const override;
    void press(Vec2 p) override;
    void release(Vec2 p) override;
    bool cancel() override;

    static int node(int column, int row) { return column * kRows + row; }
    static Vec2 node_position(int n);
    // Points fixed in rock (the left cliff's corner and face, the right cliff's face): a beam end
    // there is pinned to the world.
    static bool is_anchor(int n);
    // Adjacent points only: next to each other along the grid or across one diagonal.
    static bool can_connect(int a, int b);

    // Adds a beam between two points. Refused during a run, past the beam limit, for points that are
    // not neighbours, or if that beam is already there.
    bool add_beam(int a, int b);
    bool remove_beam(int a, int b);
    const Body& chassis() const { return world.bodies[static_cast<size_t>(car.chassis)]; }

protected:
    void clear_setup() override { beams.clear(); }
    void build() override;
    void evaluate() override;

private:
    int node_near(Vec2 p) const;
    int beam_near(Vec2 p) const;

    int from_ = -1;           // the point a new beam is being pulled from, or -1
    size_t first_beam_ = 0;   // world.bodies index of beams[0]
};

}  // namespace puzzle
