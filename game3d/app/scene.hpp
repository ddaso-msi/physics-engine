#pragma once
// Draws a level: the ground, every body coloured by its role, joints, and the marks the level asks for.
#include "painter.hpp"

namespace puzzle3d::app {

struct RoleColors {
    SDL_Color fill, edge;
};
RoleColors role_colors(Role role);

void draw_level(SDL_Renderer* ren, const demo3d::Camera& cam, const Level& level, const Ray* pointer, bool show_contacts);

}  // namespace puzzle3d::app
