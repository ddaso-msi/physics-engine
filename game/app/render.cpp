#include "render.hpp"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>

namespace puzzle::app {

using phys::Body;
using phys::Vec2;

namespace {

constexpr float kWorldH = static_cast<float>(kHeight) / kPixelsPerMeter;
constexpr int kSegments = 32;

SDL_FColor to_fcolor(SDL_Color c) {
    return {static_cast<float>(c.r) / 255.0f, static_cast<float>(c.g) / 255.0f, static_cast<float>(c.b) / 255.0f,
            static_cast<float>(c.a) / 255.0f};
}

// Convex outline as a triangle fan around its first point.
void fill_fan(SDL_Renderer* ren, const SDL_FPoint* pts, int count, SDL_Color fill) {
    SDL_Vertex verts[kSegments];
    int indices[(kSegments - 2) * 3];
    const SDL_FColor c = to_fcolor(fill);
    for (int i = 0; i < count; ++i) verts[i] = {pts[i], c, {0, 0}};
    for (int i = 0; i < count - 2; ++i) {
        indices[i * 3] = 0;
        indices[i * 3 + 1] = i + 1;
        indices[i * 3 + 2] = i + 2;
    }
    SDL_RenderGeometry(ren, nullptr, verts, count, indices, (count - 2) * 3);
}

void circle_points(SDL_FPoint centre, float radius, SDL_FPoint* out) {
    for (int i = 0; i < kSegments; ++i) {
        const float a = 2.0f * phys::kPi * static_cast<float>(i) / kSegments;
        out[i] = {centre.x + radius * std::cos(a), centre.y + radius * std::sin(a)};
    }
}

}  // namespace

SDL_FPoint to_screen(Vec2 w) { return {w.x * kPixelsPerMeter, (kWorldH - w.y) * kPixelsPerMeter}; }
Vec2 to_world(float sx, float sy) { return {sx / kPixelsPerMeter, kWorldH - sy / kPixelsPerMeter}; }

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
    SDL_FPoint pts[kSegments];
    circle_points(centre, radius, pts);
    fill_fan(ren, pts, kSegments, fill);
}

void draw_body(SDL_Renderer* ren, const Body& b, SDL_Color fill, SDL_Color edge) {
    SDL_FPoint pts[kSegments + 1];
    int count = 0;
    if (b.shape.type == phys::Shape::Type::Circle) {
        circle_points(to_screen(b.pos), b.shape.circle.radius * kPixelsPerMeter, pts);
        count = kSegments;
    } else {
        const phys::Polygon& p = b.shape.polygon;
        for (int i = 0; i < p.count; ++i) pts[i] = to_screen(apply(b.transform(), p.vertices[i]));
        count = p.count;
    }
    fill_fan(ren, pts, count, fill);
    pts[count] = pts[0];
    set_color(ren, edge);
    SDL_RenderLines(ren, pts, count + 1);
    if (b.shape.type == phys::Shape::Type::Circle) {
        const SDL_FPoint c = to_screen(b.pos), tip = to_screen(b.pos + b.q.x_axis() * b.shape.circle.radius);
        SDL_RenderLine(ren, c.x, c.y, tip.x, tip.y);
    }
}

void dashed_line(SDL_Renderer* ren, Vec2 from, Vec2 to, SDL_Color c) {
    set_color(ren, c);
    constexpr float kDash = 0.2f;
    const float length = distance(from, to);
    const Vec2 dir = (to - from).normalized();
    for (float d = 0; d < length; d += 2 * kDash) {
        const SDL_FPoint a = to_screen(from + dir * d), b = to_screen(from + dir * std::fmin(d + kDash, length));
        SDL_RenderLine(ren, a.x, a.y, b.x, b.y);
    }
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
    text_centered(ren, b.rect.x + b.rect.w * 0.5f, b.rect.y + (b.rect.h - 16.0f) * 0.5f, 2.0f,
                  enabled ? color::kText : color::kDim, b.label);
}

}  // namespace puzzle::app
