#include "hud.hpp"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>

namespace puzzle3d::app {

void set_color(SDL_Renderer* ren, SDL_Color c) { SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, c.a); }

void fill_rect(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color fill) {
    set_color(ren, fill);
    SDL_RenderFillRect(ren, &r);
}

void outline_rect(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color edge) {
    set_color(ren, edge);
    SDL_RenderRect(ren, &r);
}

void fill_disc(SDL_Renderer* ren, SDL_FPoint centre, float radius, SDL_Color fill) {
    constexpr int kSegments = 24;
    SDL_Vertex verts[kSegments + 1];
    int indices[kSegments * 3];
    const SDL_FColor c{static_cast<float>(fill.r) / 255.0f, static_cast<float>(fill.g) / 255.0f, static_cast<float>(fill.b) / 255.0f,
                       static_cast<float>(fill.a) / 255.0f};
    verts[0] = {centre, c, {0, 0}};
    for (int i = 0; i < kSegments; ++i) {
        const float a = 6.2831853f * static_cast<float>(i) / kSegments;
        verts[i + 1] = {{centre.x + radius * std::cos(a), centre.y + radius * std::sin(a)}, c, {0, 0}};
        indices[i * 3] = 0;
        indices[i * 3 + 1] = i + 1;
        indices[i * 3 + 2] = (i + 1) % kSegments + 1;
    }
    SDL_RenderGeometry(ren, nullptr, verts, kSegments + 1, indices, kSegments * 3);
}

void text(SDL_Renderer* ren, float x, float y, float scale, SDL_Color c, const char* fmt, ...) {
    char buf[256];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    set_color(ren, c);
    SDL_SetRenderScale(ren, scale, scale);
    SDL_RenderDebugText(ren, x / scale, y / scale, buf);
    SDL_SetRenderScale(ren, 1.0f, 1.0f);
}

float text_width(const char* str, float scale) {
    return static_cast<float>(std::strlen(str)) * SDL_DEBUG_TEXT_FONT_CHARACTER_SIZE * scale;
}

void text_centered(SDL_Renderer* ren, float cx, float y, float scale, SDL_Color c, const char* str) {
    text(ren, cx - text_width(str, scale) * 0.5f, y, scale, c, "%s", str);
}

void draw_button(SDL_Renderer* ren, const Button& b, bool hovered, bool enabled) {
    fill_rect(ren, b.rect, hovered && enabled ? SDL_Color{52, 60, 82, 255} : SDL_Color{36, 41, 56, 255});
    outline_rect(ren, b.rect, hovered && enabled ? color::kToolEdge : color::kPanelEdge);
    text_centered(ren, b.rect.x + b.rect.w * 0.5f, b.rect.y + (b.rect.h - 16.0f) * 0.5f, 2.0f, enabled ? color::kText : color::kDim, b.label);
}

void draw_panel(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color edge) {
    fill_rect(ren, r, color::kPanel);
    outline_rect(ren, r, edge);
    outline_rect(ren, {r.x + 1, r.y + 1, r.w - 2, r.h - 2}, edge);
}

}  // namespace puzzle3d::app
