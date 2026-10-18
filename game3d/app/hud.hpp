#pragma once
// Flat drawing for the interface: rectangles, discs, scaled text and buttons, in window pixels.
#include <SDL3/SDL.h>

namespace puzzle3d::app {

constexpr int kWidth = 1280, kHeight = 720;

namespace color {
constexpr SDL_Color kSky{30, 35, 48, 255}, kSkyLow{48, 55, 72, 255};
constexpr SDL_Color kText{226, 230, 240, 255}, kDim{140, 148, 170, 255};
constexpr SDL_Color kPanel{18, 21, 30, 235}, kPanelEdge{70, 78, 100, 255};
constexpr SDL_Color kGood{110, 220, 150, 255}, kBad{240, 110, 100, 255};
constexpr SDL_Color kTarget{242, 168, 59, 255}, kToolEdge{170, 222, 255, 255};
}  // namespace color

void set_color(SDL_Renderer* ren, SDL_Color c);
void fill_rect(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color fill);
void outline_rect(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color edge);
void fill_disc(SDL_Renderer* ren, SDL_FPoint centre, float radius, SDL_Color fill);

// Text in SDL's built-in 8x8 font, enlarged `scale` times. (x, y) is the top-left corner in pixels.
void text(SDL_Renderer* ren, float x, float y, float scale, SDL_Color c, SDL_PRINTF_FORMAT_STRING const char* fmt, ...)
    SDL_PRINTF_VARARG_FUNC(6);
void text_centered(SDL_Renderer* ren, float cx, float y, float scale, SDL_Color c, const char* str);
float text_width(const char* str, float scale);

struct Button {
    SDL_FRect rect;
    const char* label;
    bool hit(float x, float y) const { return x >= rect.x && x < rect.x + rect.w && y >= rect.y && y < rect.y + rect.h; }
};
void draw_button(SDL_Renderer* ren, const Button& b, bool hovered, bool enabled = true);
void draw_panel(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color edge);

}  // namespace puzzle3d::app
