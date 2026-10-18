#include "scene.hpp"

namespace puzzle3d::app {

RoleColors role_colors(Role role) {
    switch (role) {
        case Role::Ground: return {{70, 76, 92, 255}, {120, 128, 150, 255}};
        case Role::Scenery: return {{112, 120, 144, 255}, {170, 178, 200, 255}};
        case Role::Target: return {{242, 160, 52, 255}, {255, 214, 140, 255}};
        case Role::Done: return {{80, 170, 120, 255}, {130, 220, 170, 255}};
        case Role::Tool: return {{70, 160, 236, 255}, {170, 222, 255, 255}};
        case Role::Load: return {{168, 128, 236, 255}, {222, 204, 255, 255}};
        case Role::Selected: return {{140, 214, 255, 255}, {255, 255, 255, 255}};
        case Role::Invalid: return {{214, 80, 72, 255}, {255, 140, 130, 255}};
        case Role::Guide: return {{120, 128, 150, 255}, {200, 206, 222, 255}};
    }
    return {{112, 120, 144, 255}, {170, 178, 200, 255}};
}

void draw_level(SDL_Renderer* ren, const demo3d::Camera& cam, const Level& level, const Ray* pointer, bool show_contacts) {
    Painter paint(ren, cam);
    const World& world = level.world;
    Overlay marks;
    level.overlay(marks, pointer);

    paint.sky();
    if (!world.bodies.empty()) paint.ground(world.bodies[0]);
    for (size_t i = 1; i < world.bodies.size(); ++i)
        if (world.bodies[i].type == phys3d::BodyType::Dynamic) paint.shadow(world.bodies[i]);
    for (const Overlay::Ghost& g : marks.ghosts) paint.shadow(g.body);

    for (const Overlay::Zone& z : marks.zones) paint.zone_floor(z.lo, z.hi, role_colors(z.role).fill);
    for (size_t i = 1; i < world.bodies.size(); ++i) paint.body(world.bodies[i], {role_colors(level.role(i)).fill, 255});
    for (const Overlay::Ghost& g : marks.ghosts) paint.body(g.body, {role_colors(g.role).fill, 120});
    paint.flush();

    for (const Overlay::Zone& z : marks.zones) {  // the zone's edges
        const SDL_Color c = role_colors(z.role).edge;
        const Vec3 lo = z.lo, hi = z.hi;
        const Vec3 b[4] = {{lo.x, lo.y, lo.z}, {hi.x, lo.y, lo.z}, {hi.x, lo.y, hi.z}, {lo.x, lo.y, hi.z}};
        for (int i = 0; i < 4; ++i) {
            const Vec3 top_a{b[i].x, hi.y, b[i].z}, top_b{b[(i + 1) % 4].x, hi.y, b[(i + 1) % 4].z};
            paint.line(b[i], b[(i + 1) % 4], c);
            paint.line(top_a, top_b, c, true, 1.0f);
            paint.line(b[i], top_a, c, true, 1.0f);
        }
    }
    for (const phys3d::Joint& j : world.joints) {
        const Vec3 pa = j.world_anchor_a(world.bodies), pb = j.world_anchor_b(world.bodies);
        const SDL_Color c{236, 240, 248, 255};
        if (j.type == phys3d::JointType::Distance) paint.line(pa, pb, c);
        if (j.a < 0 || j.b < 0) paint.marker(j.a < 0 ? pa : pb, 0.07f, c);
        else paint.marker(pa, 0.05f, {236, 240, 248, 200});
    }
    for (const Overlay::Line& l : marks.lines) paint.line(l.a, l.b, role_colors(l.role).edge, l.dashed);
    for (const Overlay::Dot& d : marks.dots) paint.marker(d.p, d.radius, role_colors(d.role).edge);
    for (const Overlay::Label& l : marks.labels) paint.label(l.p, l.text.c_str(), role_colors(l.role).edge);

    if (show_contacts)
        for (const phys3d::ContactPair& c : world.contacts())
            for (int k = 0; k < c.manifold.count; ++k) {
                const Vec3 p = c.manifold.points[k].point;
                paint.line(p, p + c.manifold.normal * 0.3f, color::kGood);
                paint.marker(p, 0.03f, {255, 240, 150, 255});
            }
}

}  // namespace puzzle3d::app
