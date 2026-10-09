#include "scene.hpp"

namespace puzzle::app {

RoleColors role_colors(Role role) {
    switch (role) {
        case Role::Scenery: return {color::kStatic, color::kStaticEdge};
        case Role::Target: return {color::kTarget, color::kTargetEdge};
        case Role::Done: return {color::kDone, color::kDoneEdge};
        case Role::Tool: return {color::kTool, color::kToolEdge};
        case Role::Load: return {{168, 132, 232, 255}, {222, 204, 255, 255}};
        case Role::Selected: return {{120, 200, 250, 255}, {255, 255, 255, 255}};
        case Role::Invalid: return {{150, 60, 60, 255}, color::kBad};
        case Role::Guide: return {{96, 104, 128, 255}, {190, 196, 214, 255}};
    }
    return {color::kStatic, color::kStaticEdge};
}

void draw_grid(SDL_Renderer* ren) {
    set_color(ren, color::kGrid);
    for (int x = 1; x < kWidth / static_cast<int>(kPixelsPerMeter); ++x) {
        const float sx = static_cast<float>(x) * kPixelsPerMeter;
        SDL_RenderLine(ren, sx, 0, sx, static_cast<float>(kHeight));
    }
    for (int y = 1; y <= kHeight / static_cast<int>(kPixelsPerMeter); ++y) {
        const SDL_FPoint p = to_screen({0, static_cast<float>(y)});
        SDL_RenderLine(ren, 0, p.y, static_cast<float>(kWidth), p.y);
    }
}

namespace {

void draw_joint(SDL_Renderer* ren, const World& world, const phys::Joint& j) {
    const SDL_FPoint a = to_screen(j.world_anchor_a(world.bodies)), b = to_screen(j.world_anchor_b(world.bodies));
    if (j.type == phys::JointType::Distance) {  // a rod or rope, with its fixed end marked
        set_color(ren, color::kText);
        SDL_RenderLine(ren, a.x, a.y, b.x, b.y);
        if (j.a < 0) fill_disc(ren, a, 5.0f, color::kText);
        if (j.b < 0) fill_disc(ren, b, 5.0f, color::kText);
    } else if (j.type == phys::JointType::Revolute) {  // a pin: solid when it is fixed to the world
        fill_disc(ren, a, 5.0f, j.a < 0 || j.b < 0 ? color::kText : color::kBackground);
        fill_disc(ren, a, 2.5f, j.a < 0 || j.b < 0 ? color::kBackground : color::kText);
    }
}

SDL_Color faded(SDL_Color c, Uint8 alpha) { return {c.r, c.g, c.b, alpha}; }

}  // namespace

void draw_level(SDL_Renderer* ren, const Level& level, Vec2 pointer, bool pointer_active, bool show_contacts) {
    draw_grid(ren);
    Overlay marks;
    level.overlay(marks, pointer, pointer_active);

    for (const Overlay::Zone& z : marks.zones) {  // under the bodies, so what is inside a zone stays readable
        const SDL_FPoint lo = to_screen(z.lo), hi = to_screen(z.hi);
        const SDL_FRect r{lo.x, hi.y, hi.x - lo.x, lo.y - hi.y};
        fill_rect(ren, r, faded(role_colors(z.role).fill, 38));
        outline_rect(ren, r, faded(role_colors(z.role).edge, 150));
    }

    const World& world = level.world;
    for (size_t i = 0; i < world.bodies.size(); ++i) {
        const RoleColors c = role_colors(level.role(i));
        draw_body(ren, world.bodies[i], c.fill, c.edge);
    }
    for (const phys::Joint& j : world.joints) draw_joint(ren, world, j);

    for (const Overlay::Ghost& g : marks.ghosts) {
        const RoleColors c = role_colors(g.role);
        draw_body(ren, g.body, faded(c.fill, 110), c.edge);
    }
    for (const Overlay::Line& l : marks.lines) {
        const SDL_Color c = role_colors(l.role).edge;
        if (l.dashed) {
            dashed_line(ren, l.a, l.b, c);
        } else {
            const SDL_FPoint a = to_screen(l.a), b = to_screen(l.b);
            set_color(ren, c);
            SDL_RenderLine(ren, a.x, a.y, b.x, b.y);
            SDL_RenderLine(ren, a.x, a.y + 1, b.x, b.y + 1);
        }
    }
    for (const Overlay::Dot& d : marks.dots) {
        const RoleColors c = role_colors(d.role);
        const float r = d.radius * kPixelsPerMeter;
        if (r > 6.0f) fill_disc(ren, to_screen(d.p), r + 2.0f, c.edge);
        fill_disc(ren, to_screen(d.p), r, r > 6.0f ? c.fill : c.edge);
    }
    for (const Overlay::Label& l : marks.labels) {
        const SDL_FPoint p = to_screen(l.p);
        text(ren, p.x, p.y, 2.0f, role_colors(l.role).edge, "%s", l.text.c_str());
    }

    if (show_contacts)
        for (const phys::ContactPair& c : world.contacts())
            for (int k = 0; k < c.manifold.count; ++k) {
                const SDL_FPoint p = to_screen(c.manifold.points[k].point);
                const SDL_FPoint tip = to_screen(c.manifold.points[k].point + c.manifold.normal * 0.3f);
                set_color(ren, color::kGood);
                SDL_RenderLine(ren, p.x, p.y, tip.x, tip.y);
                fill_rect(ren, {p.x - 2.0f, p.y - 2.0f, 4.0f, 4.0f}, {255, 240, 150, 255});
            }
}

}  // namespace puzzle::app
