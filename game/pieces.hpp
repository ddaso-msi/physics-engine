#pragma once
// Levels where the player places loose pieces during setup. A piece starts in a tray, is dragged into
// the build area, and drops straight down until it rests on something. While the world is frozen
// nothing can push anything, so a piece is only accepted where it overlaps nothing.
#include "level.hpp"

#include <phys/aabb.hpp>

#include <vector>

namespace puzzle {

struct PieceKind {
    const char* label;
    Vec2 half;        // half extents of a box...
    Real radius = 0;  // ...or, if positive, a disc of this radius
    Real density = 1, friction = 0.6f, restitution = 0.05f;
};

struct Piece {
    PieceKind kind;
    Vec2 home;          // where it waits in the tray
    Vec2 pos;
    Real angle = 0;
    bool placed = false;
    int body = -1;      // index into world.bodies while placed, else -1
};

class PieceLevel : public Level {
public:
    std::vector<Piece> pieces;

    Role role(size_t body) const override;
    void overlay(Overlay& out, Vec2 pointer, bool pointer_active) const override;
    void press(Vec2 p) override;
    void drag(Vec2 p) override;
    void release(Vec2 p) override;
    void rotate(int direction) override;
    bool cancel() override;

    // Puts a piece at `pos` (it then drops until it rests on something). Refused, leaving everything as
    // it was, when the level is running or the spot is outside the build area or overlaps something.
    bool place(int piece, Vec2 pos, Real angle = 0);
    // Sends a placed piece back to the tray.
    void unplace(int piece);
    int placed_count() const;
    std::vector<int> placed_bodies() const;
    bool holding() const { return held_ >= 0; }

protected:
    virtual std::vector<Piece> initial_pieces() const = 0;
    virtual void build_scenery() = 0;       // everything that is not a player's piece
    virtual phys::AABB build_area() const = 0;
    virtual Real rotation_step() const { return phys::kPi / 12; }
    virtual Role scenery_role(size_t body) const;

    void clear_setup() override;
    void build() override;

private:
    bool fits(int piece) const;  // inside the build area and overlapping nothing
    void settle(int piece);      // lower it until it touches something
    int piece_at(Vec2 p) const;

    int held_ = -1;
    bool held_fits_ = false;
    Vec2 grab_;      // from the pointer to the held piece's centre
    Piece before_;   // the held piece as it was, to put it back
};

}  // namespace puzzle
