#pragma once
// Draws a level: its bodies coloured by role, its joints, and the marks it asks for.
#include "render.hpp"

#include <level.hpp>

namespace puzzle::app {

struct RoleColors {
    SDL_Color fill, edge;
};
RoleColors role_colors(Role role);

void draw_grid(SDL_Renderer* ren);
void draw_level(SDL_Renderer* ren, const Level& level, Vec2 pointer, bool pointer_active, bool show_contacts);

}  // namespace puzzle::app
