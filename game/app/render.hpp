#pragma once
// Drawing helpers for the game: world-to-screen mapping, filled shapes, scaled text, buttons.
#include <SDL3/SDL.h>
#include <phys/world.hpp>

namespace puzzle::app {

constexpr int kWidth = 1280, kHeight = 720;
constexpr float kPixelsPerMeter = 64.0f;  // the window shows the 20 x 11.25 m arena exactly

SDL_FPoint to_screen(phys::Vec2 w);
phys::Vec2 to_world(float sx, float sy);

// The palette. Targets are the only warm colour on screen.
namespace color {
constexpr SDL_Color kBackground{27, 31, 42, 255}, kGrid{34, 39, 52, 255};
constexpr SDL_Color kStatic{58, 65, 84, 255}, kStaticEdge{126, 136, 162, 255};
constexpr SDL_Color kTarget{242, 168, 59, 255}, kTargetEdge{255, 214, 140, 255};
constexpr SDL_Color kDone{70, 122, 100, 255}, kDoneEdge{120, 200, 160, 255};
constexpr SDL_Color kTool{74, 170, 232, 255}, kToolEdge{170, 222, 255, 255};
constexpr SDL_Color kText{226, 230, 240, 255}, kDim{140, 148, 170, 255};
constexpr SDL_Color kPanel{18, 21, 30, 235}, kPanelEdge{70, 78, 100, 255};
constexpr SDL_Color kGood{110, 220, 150, 255}, kBad{240, 110, 100, 255};
}  // namespace color

void set_color(SDL_Renderer* ren, SDL_Color c);
void fill_rect(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color fill);
void outline_rect(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color edge);
// A disc given in screen pixels.
void fill_disc(SDL_Renderer* ren, SDL_FPoint centre, float radius, SDL_Color fill);
// A body, filled and outlined. Circles get a spoke so spin is visible.
void draw_body(SDL_Renderer* ren, const phys::Body& b, SDL_Color fill, SDL_Color edge);
// A dashed line between two world points.
void dashed_line(SDL_Renderer* ren, phys::Vec2 from, phys::Vec2 to, SDL_Color c);

// Text in SDL's built-in 8x8 font, enlarged `scale` times. (x, y) is the top-left corner in pixels.
void text(SDL_Renderer* ren, float x, float y, float scale, SDL_Color c, SDL_PRINTF_FORMAT_STRING const char* fmt, ...)
    SDL_PRINTF_VARARG_FUNC(6);
void text_centered(SDL_Renderer* ren, float cx, float y, float scale, SDL_Color c, const char* str);
float text_width(const char* str, float scale);

struct Button {
    SDL_FRect rect;
    const char* label;
    bool hit(float x, float y) const {
        return x >= rect.x && x < rect.x + rect.w && y >= rect.y && y < rect.y + rect.h;
    }
};
void draw_button(SDL_Renderer* ren, const Button& b, bool hovered, bool enabled = true);

}  // namespace puzzle::app
