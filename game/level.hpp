#pragma once
// Shared by every level: the arena, the fixed step, and what a level reports about its outcome.
#include <phys/world.hpp>

namespace puzzle {

using phys::Body;
using phys::Real;
using phys::Vec2;
using phys::World;

// The arena, in metres. The app shows exactly this area, so levels can be laid out in world units.
constexpr Real kWorldW = 20.0f, kWorldH = 11.25f;
constexpr Real kFloorTop = 0.5f;
// Every level advances by this step and nothing else, so a replay of the same inputs is identical.
constexpr Real kTimeStep = 1.0f / 120.0f;

enum class Status { Playing, Won, Failed };

// Floor (top at kFloorTop), side walls, and a ceiling just out of view.
void add_room(World& w);

}  // namespace puzzle
