#pragma once
// Levels where the player places loose pieces during setup. A piece waits beside the build area, is
// picked with the mouse, and while held it hangs over whatever the pointer is over, already dropped as far
// as it would fall: let go and that is where it stays. While the world is frozen nothing can push anything,
// so a piece is only accepted where it overlaps nothing and lies inside the build area.
#include "level.hpp"

#include <phys3d/aabb.hpp>

#include <vector>

namespace puzzle3d {

struct PieceKind {
    const char* label;
    Vec3 half;        // half extents of a box...
    Real radius = 0;  // ...or, if positive, a ball of this radius
    Real density = 1, friction = 0.6f, restitution = 0.05f;
};

struct Piece {
    PieceKind kind;
    Vec3 home;           // where it waits, beside the build area
    Vec3 pos;
    Real yaw = 0;        // turned about the vertical
    bool upright = false;  // stood on end: its own x axis made vertical
    bool placed = false;
    int body = -1;       // index into world.bodies while placed, else -1

    Quat orientation() const;
};

Body make_body(const Piece& p, Vec3 pos);

class PieceLevel : public Level {
public:
    std::vector<Piece> pieces;

    Role role(size_t body) const override;
    void overlay(Overlay& out, const Ray* pointer) const override;
    void press(const Ray& ray) override;
    void drag(const Ray& ray) override;
    void release(const Ray& ray) override;
    void rotate(int direction) override;
    void tip() override;
    bool cancel() override;

    // Drops a piece at (x, z) with the given turn: it falls from high up until it rests on something.
    // Refused, leaving everything as it was, during a run or where it would end outside the build area.
    bool place(int piece, Real x, Real z, Real yaw = 0, bool upright = false);
    void unplace(int piece);  // back beside the build area
    int placed_count() const;
    std::vector<int> placed_bodies() const;
    bool holding() const { return held_ >= 0; }

protected:
    virtual std::vector<Piece> initial_pieces() const = 0;
    virtual void build_scenery() = 0;  // everything that is not a player's piece
    virtual phys3d::AABB build_area() const = 0;
    virtual Real rotation_step() const { return phys3d::Real(3.14159265f / 12); }
    virtual Role scenery_role(size_t body) const;

    void clear_setup() override;
    void build() override;

private:
    bool fits(int piece) const;          // inside the build area and overlapping nothing
    void drop(int piece, Real x, Real z);  // lower it from high up until it is about to touch something
    int piece_under(const Ray& ray) const;

    int held_ = -1;
    bool held_fits_ = false;
    Piece before_;  // the held piece as it was, to put it back
    Vec3 last_spot_;
};

}  // namespace puzzle3d
