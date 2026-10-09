#include "pieces.hpp"

#include <phys/collision.hpp>

#include <algorithm>

namespace puzzle {

namespace {

Body make_body(const Piece& p, Vec2 pos) {
    const phys::Shape shape = p.kind.radius > 0 ? phys::Shape::make_circle(p.kind.radius)
                                                : phys::Shape::make_polygon(phys::Polygon::box(p.kind.half.x, p.kind.half.y));
    Body b(shape, pos, p.angle, phys::BodyType::Dynamic, p.kind.density);
    b.friction = p.kind.friction;
    b.restitution = p.kind.restitution;
    return b;
}

// Does body `index` overlap any other body by more than `tolerance`?
bool overlaps(const World& w, int index, Real tolerance) {
    const Body& body = w.bodies[static_cast<size_t>(index)];
    for (size_t i = 0; i < w.bodies.size(); ++i) {
        if (static_cast<int>(i) == index || !should_collide(body, w.bodies[i])) continue;
        phys::Manifold m;
        if (collide(body, w.bodies[i], m) && m.depth > tolerance) return true;
    }
    return false;
}

}  // namespace

void PieceLevel::clear_setup() {
    pieces = initial_pieces();
    for (Piece& p : pieces) p.pos = p.home;
    held_ = -1;
}

void PieceLevel::build() {
    build_scenery();
    for (Piece& p : pieces) p.body = p.placed ? world.add(make_body(p, p.pos)) : -1;
}

Role PieceLevel::scenery_role(size_t body) const {
    return world.bodies[body].type == phys::BodyType::Static ? Role::Scenery : Role::Guide;
}

Role PieceLevel::role(size_t body) const {
    for (size_t i = 0; i < pieces.size(); ++i)
        if (pieces[i].body == static_cast<int>(body)) {
            if (static_cast<int>(i) != held_) return Role::Tool;
            return held_fits_ ? Role::Selected : Role::Invalid;
        }
    return scenery_role(body);
}

void PieceLevel::overlay(Overlay& out, Vec2, bool) const {
    if (pieces.empty()) return;
    phys::AABB tray{pieces[0].home, pieces[0].home};
    for (const Piece& p : pieces) {
        Piece at_rest = p;
        at_rest.angle = 0;
        const Body ghost = make_body(at_rest, p.home);
        tray = merge(tray, compute_aabb(ghost));
        if (!p.placed) out.ghosts.push_back({ghost, phase == Phase::Setup ? Role::Tool : Role::Guide});
    }
    tray = fatten(tray, 0.3f);
    out.zones.push_back({tray.lo, tray.hi, Role::Guide});
    out.labels.push_back({{tray.lo.x, tray.hi.y + 0.4f}, "pieces", Role::Guide});
    if (phase == Phase::Setup) {
        const phys::AABB area = build_area();
        out.zones.push_back({area.lo, area.hi, Role::Tool});
    }
}

int PieceLevel::placed_count() const {
    return static_cast<int>(std::count_if(pieces.begin(), pieces.end(), [](const Piece& p) { return p.placed; }));
}

std::vector<int> PieceLevel::placed_bodies() const {
    std::vector<int> out;
    for (const Piece& p : pieces)
        if (p.body >= 0) out.push_back(p.body);
    return out;
}

bool PieceLevel::fits(int piece) const {
    const Piece& p = pieces[static_cast<size_t>(piece)];
    if (p.body < 0) return false;
    if (!contains(build_area(), compute_aabb(world.bodies[static_cast<size_t>(p.body)]))) return false;
    return !overlaps(world, p.body, 0.01f);
}

void PieceLevel::settle(int piece) {
    Piece& p = pieces[static_cast<size_t>(piece)];
    Body& body = world.bodies[static_cast<size_t>(p.body)];
    // Coarse steps down to the first touch, then back up and repeat with fine ones: it ends within a
    // couple of millimetres of whatever is below, touching nothing.
    for (const Real step : {0.02f, 0.002f}) {
        while (!overlaps(world, p.body, 0) && body.pos.y > 0) body.pos.y -= step;
        body.pos.y += step;
    }
    p.pos = body.pos;
}

bool PieceLevel::place(int piece, Vec2 pos, Real angle) {
    if (phase != Phase::Setup || piece < 0 || piece >= static_cast<int>(pieces.size())) return false;
    cancel();
    Piece& p = pieces[static_cast<size_t>(piece)];
    const Piece before = p;
    p.pos = pos;
    p.angle = angle;
    p.placed = true;
    rebuild();
    if (!fits(piece)) {
        p = before;
        rebuild();
        return false;
    }
    settle(piece);
    return true;
}

void PieceLevel::unplace(int piece) {
    if (phase != Phase::Setup) return;
    cancel();
    Piece& p = pieces[static_cast<size_t>(piece)];
    p.placed = false;
    p.pos = p.home;
    p.angle = 0;
    rebuild();
}

int PieceLevel::piece_at(Vec2 point) const {
    for (int i = static_cast<int>(pieces.size()) - 1; i >= 0; --i) {
        const Piece& p = pieces[static_cast<size_t>(i)];
        if (p.placed ? world.bodies[static_cast<size_t>(p.body)].contains(point) : make_body(p, p.home).contains(point)) return i;
    }
    return -1;
}

void PieceLevel::press(Vec2 point) {
    if (phase != Phase::Setup || held_ >= 0) return;
    const int i = piece_at(point);
    if (i < 0) return;
    held_ = i;
    Piece& p = pieces[static_cast<size_t>(i)];
    before_ = p;
    grab_ = p.pos - point;
    p.placed = true;
    drag(point);
}

void PieceLevel::drag(Vec2 point) {
    if (held_ < 0) return;
    pieces[static_cast<size_t>(held_)].pos = point + grab_;
    rebuild();
    held_fits_ = fits(held_);
}

void PieceLevel::release(Vec2 point) {
    if (held_ < 0) return;
    drag(point);
    const int i = held_;
    if (!held_fits_) {  // dropped somewhere it cannot go: back where it came from
        cancel();
        return;
    }
    held_ = -1;
    settle(i);
}

void PieceLevel::rotate(int direction) {
    if (held_ < 0) return;
    Piece& p = pieces[static_cast<size_t>(held_)];
    p.angle += static_cast<Real>(direction) * rotation_step();
    rebuild();
    held_fits_ = fits(held_);
}

bool PieceLevel::cancel() {
    if (held_ < 0) return false;
    pieces[static_cast<size_t>(held_)] = before_;
    held_ = -1;
    rebuild();
    return true;
}

}  // namespace puzzle
