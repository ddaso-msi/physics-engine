#pragma once
// What every level has in common: the arena, the fixed step, the setup/run life cycle, and the plain
// data a level hands to whoever draws it. Nothing here knows about windows or input devices.
#include <phys/world.hpp>

#include <span>
#include <string>
#include <vector>

namespace puzzle {

using phys::Body;
using phys::Real;
using phys::Vec2;
using phys::World;

// The arena, in metres. The app shows exactly this area, so levels are laid out in world units.
constexpr Real kWorldW = 20.0f, kWorldH = 11.25f;
constexpr Real kFloorTop = 0.5f;
// Every level advances by this step and nothing else, so a replay of the same inputs is identical.
constexpr Real kTimeStep = 1.0f / 120.0f;

enum class Status { Playing, Won, Failed };
// Setup: the world is built but frozen, and the player arranges things. Running: physics advances.
enum class Phase { Setup, Running };

// What a body (or a mark on the screen) means to the player. The app picks colours from this.
enum class Role { Scenery, Target, Done, Tool, Load, Selected, Invalid, Guide };

// Marks a level wants drawn on top of the world, in world coordinates.
struct Overlay {
    struct Line { Vec2 a, b; Role role; bool dashed; };
    struct Zone { Vec2 lo, hi; Role role; };
    struct Dot { Vec2 p; Real radius; Role role; };
    struct Label { Vec2 p; std::string text; Role role; };
    struct Ghost { Body body; Role role; };  // drawn like a body but not part of the world
    std::vector<Line> lines;
    std::vector<Zone> zones;
    std::vector<Dot> dots;
    std::vector<Label> labels;
    std::vector<Ghost> ghosts;
};

// Floor (top at kFloorTop), side walls, and a ceiling just out of view.
void add_room(World& w);

class Level {
public:
    World world;
    Phase phase = Phase::Setup;
    Status status = Status::Playing;
    double sim_time = 0;  // simulated seconds in the current run
    double win_time = 0;  // sim_time at which the objective was met
    int attempt = 1;      // counts runs that were actually tried; survives reset()

    virtual ~Level() = default;

    virtual const char* name() const = 0;
    virtual const char* objective() const = 0;
    virtual std::span<const char* const> help() const = 0;  // how to play, a line each
    virtual const char* hint() const = 0;                   // what to do right now
    virtual std::string stats() const = 0;                  // counters for the status bar
    virtual std::string result() const = 0;                 // the score line shown on a win
    // Levels without a setup phase run all the time (the player acts while physics is live).
    virtual bool has_setup() const { return true; }
    virtual Role role(size_t body) const;
    virtual void overlay(Overlay& out, Vec2 pointer, bool pointer_active) const;

    // Pointer input in world coordinates. press/release come once per click; drag only between them.
    virtual void press(Vec2) {}
    virtual void drag(Vec2) {}
    virtual void release(Vec2) {}
    virtual void rotate(int /*direction*/) {}
    // Abandons an interaction in progress. Returns false if there was nothing to cancel.
    virtual bool cancel() { return false; }

    // The level exactly as it first appeared: the player's setup is forgotten too.
    void reset();
    // Setup -> Running, from the player's current setup.
    void start();
    // Running -> Setup, keeping the player's setup so it can be revised.
    void stop();
    // One fixed step of physics and then the objective. Does nothing during setup.
    void step(Real dt = kTimeStep);

protected:
    virtual void clear_setup() = 0;  // forget what the player arranged
    virtual void build() = 0;        // fill the (empty) world from the level and the player's setup
    virtual void evaluate() = 0;     // called after every step while Playing or Failed
    // A new World, then build(). Never patched in place: body indices are identities in this engine and
    // there is no removing one, so rebuilding is both the simplest and the safest way to change a scene.
    void rebuild();
    void win();
    void fail();
    void mark_tried() { tried_ = true; }

private:
    void begin_run();
    bool tried_ = false;
};

}  // namespace puzzle
