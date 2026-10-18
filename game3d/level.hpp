#pragma once
// What every 3D level has in common: the ground, the fixed step, the setup/run life cycle, and the plain
// data a level hands to whoever draws it. Nothing here knows about windows or input devices.
#include "raycast.hpp"

#include <phys3d/world.hpp>

#include <span>
#include <string>
#include <vector>

namespace puzzle3d {

using phys3d::Body;
using phys3d::Quat;
using phys3d::Real;
using phys3d::Vec3;
using phys3d::World;

// Every level advances by this step and nothing else, so a replay of the same inputs is identical.
constexpr Real kTimeStep = 1.0f / 120.0f;
constexpr Real kGroundY = 0;  // the top of the ground; y is up

enum class Status { Playing, Won, Failed };
// Setup: the world is built but frozen, and the player arranges things. Running: physics advances.
enum class Phase { Setup, Running };

// What a body (or a mark) means to the player. The app picks colours from this.
enum class Role { Ground, Scenery, Target, Done, Tool, Load, Selected, Invalid, Guide };

// Marks a level wants drawn with the world, in world coordinates.
struct Overlay {
    struct Line { Vec3 a, b; Role role; bool dashed; };
    struct Zone { Vec3 lo, hi; Role role; };  // an axis-aligned box: its floor tinted, its edges drawn
    struct Dot { Vec3 p; Real radius; Role role; };
    struct Label { Vec3 p; std::string text; Role role; };
    struct Ghost { Body body; Role role; };    // drawn like a body, see-through, not part of the world
    std::vector<Line> lines;
    std::vector<Zone> zones;
    std::vector<Dot> dots;
    std::vector<Label> labels;
    std::vector<Ghost> ghosts;
};

// Where the camera starts for a level (see demo3d::Camera for what the numbers mean).
struct View {
    Vec3 target;
    Real yaw, pitch, distance;
};

// The ground: a slab 60 m across whose top is at kGroundY. Always body 0.
void add_ground(World& w);

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
    virtual View view() const = 0;
    // Levels without a setup phase run all the time (the player acts while physics is live).
    virtual bool has_setup() const { return true; }
    virtual Role role(size_t body) const;
    // `pointer` is the ray under the mouse, or null when the mouse is not over the scene.
    virtual void overlay(Overlay& out, const Ray* pointer) const;

    // Pointer input, as a ray from the camera through the mouse. press/release come once per click;
    // drag only between them.
    virtual void press(const Ray&) {}
    virtual void drag(const Ray&) {}
    virtual void release(const Ray&) {}
    virtual void rotate(int /*direction*/) {}  // turn a held piece about the vertical
    virtual void tip() {}                      // stand a held piece on end, or lay it back down
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

}  // namespace puzzle3d
