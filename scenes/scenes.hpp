#pragma once
// Stage 9: ready-made scenes for the demo gallery and the tests. Each builder adds bodies and joints to a
// World and returns the indices a caller needs to poke at them.
#include <phys/world.hpp>

#include <vector>

namespace phys::scenes {

// A static box, handy for floors, walls and ramps.
int add_static_box(World& w, Vec2 half_extents, Vec2 pos, Real angle = 0);

// Newton's cradle: `count` equal balls hanging on ropes, just touching. Perfectly elastic, no friction.
struct Cradle {
    std::vector<int> balls;
    std::vector<Vec2> pivots;
    Real rope = 0;
};
Cradle build_cradle(World& w, Vec2 top_center, int count = 5, Real rope = 4, Real ball_radius = 0.5f);
// Swing ball `which` out to `angle` radians (positive swings right) on its rope, at rest.
void lift_cradle_ball(World& w, const Cradle& c, int which, Real angle);

// A rag doll: ten bodies and nine hinges, each with a joint limit, standing with its torso centred at
// `torso_center`. Its parts do not collide with each other (they share a collision group) but do with
// everything else. `group` (1 to 14) must differ between rag dolls and cars if they should collide.
struct Ragdoll {
    int torso = 0, head = 0;
    int upper_arm[2] = {}, forearm[2] = {}, thigh[2] = {}, shin[2] = {};
    std::vector<int> hinges;  // indices into World::joints
    int knee_hinge[2] = {};
};
Ragdoll build_ragdoll(World& w, Vec2 torso_center, int group);
// Height of the torso centre above the soles of the feet when standing.
constexpr Real kRagdollTorsoHeight = static_cast<Real>(1.74);

// A car: a chassis on two wheels, each wheel hung from the chassis by a vertical slide with end stops
// (a prismatic joint) and a spring-damper (a soft distance joint), with a motor in each wheel's hinge.
struct Car {
    int chassis = 0;
    int wheel[2] = {};         // [0] rear, [1] front
    int carrier[2] = {};       // the small massless-ish body that joins slide, spring and wheel
    int drive_hinge[2] = {};   // indices into World::joints
};
Car build_car(World& w, Vec2 chassis_center, int group);
// Positive `wheel_speed` (rad/s) drives the car to the right. max_torque 0 lets the wheels coast.
void set_car_throttle(World& w, const Car& car, Real wheel_speed, Real max_torque);
constexpr Real kCarWheelRadius = static_cast<Real>(0.45);

}  // namespace phys::scenes
