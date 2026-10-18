#include "pieces.hpp"

#include <phys3d/collision.hpp>

#include <algorithm>
#include <cmath>

namespace puzzle3d {

namespace {

constexpr Real kDropHeight = 14.0f;

// Does body `index` overlap any other body by more than `tolerance`?
bool overlaps(const World& w, int index, Real tolerance) {
    const Body& body = w.bodies[static_cast<size_t>(index)];
    for (size_t i = 0; i < w.bodies.size(); ++i) {
        if (static_cast<int>(i) == index) continue;
        phys3d::Manifold m;
        if (collide(body, w.bodies[i], m) && m.depth > tolerance) return true;
    }
    return false;
}

}  // namespace

Quat Piece::orientation() const {
    const Quat turn = Quat::from_axis_angle({0, 1, 0}, yaw);
    return upright ? turn * Quat::from_axis_angle({0, 0, 1}, 3.14159265f / 2) : turn;
}

Body make_body(const Piece& p, Vec3 pos) {
    Body b = p.kind.radius > 0 ? Body::solid_sphere(p.kind.radius, p.kind.density, pos)
                               : Body::solid_box(p.kind.half, p.kind.density, pos, p.orientation());
    b.friction = p.kind.friction;
    b.restitution = p.kind.restitution;
    return b;
}

void PieceLevel::clear_setup() {
    pieces = initial_pieces();
    for (Piece& p : pieces) p.pos = p.home;
    held_ = -1;
}

void PieceLevel::build() {
    build_scenery();
    for (Piece& p : pieces) p.body = p.placed ? world.add(make_body(p, p.pos)) : -1;
}

Role PieceLevel::scenery_role(size_t body) const { return Level::role(body); }

Role PieceLevel::role(size_t body) const {
    for (size_t i = 0; i < pieces.size(); ++i)
        if (pieces[i].body == static_cast<int>(body)) {
            if (static_cast<int>(i) != held_) return Role::Tool;
            return held_fits_ ? Role::Selected : Role::Invalid;
        }
    return scenery_role(body);
}

void PieceLevel::overlay(Overlay& out, const Ray*) const {
    for (const Piece& p : pieces)
        if (!p.placed) out.ghosts.push_back({make_body(p, p.home), phase == Phase::Setup ? Role::Tool : Role::Guide});
    if (phase == Phase::Setup) {
        const phys3d::AABB area = build_area();
        out.zones.push_back({area.lo, {area.hi.x, area.lo.y + 0.05f, area.hi.z}, Role::Tool});
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
    if (!phys3d::contains(build_area(), phys3d::compute_aabb(world.bodies[static_cast<size_t>(p.body)]))) return false;
    // drop() always stops a piece just short of touching anything, so this cannot fail today; it stays as a
    // guard for any future way of placing a piece.
    return !overlaps(world, p.body, 0.01f);
}

void PieceLevel::drop(int piece, Real x, Real z) {
    Piece& p = pieces[static_cast<size_t>(piece)];
    p.pos = {x, kDropHeight, z};
    p.placed = true;
    rebuild();
    Body& body = world.bodies[static_cast<size_t>(p.body)];
    // Coarse steps down to the first touch, then back up and repeat with fine ones: it ends within a few
    // millimetres of whatever is below, touching nothing.
    for (const Real step : {0.05f, 0.004f}) {
        while (!overlaps(world, p.body, 0) && body.pos.y > kGroundY - 1) body.pos.y -= step;
        body.pos.y += step;
    }
    p.pos = body.pos;
}

bool PieceLevel::place(int piece, Real x, Real z, Real yaw, bool upright) {
    if (phase != Phase::Setup || piece < 0 || piece >= static_cast<int>(pieces.size())) return false;
    cancel();
    Piece& p = pieces[static_cast<size_t>(piece)];
    const Piece before = p;
    p.yaw = yaw;
    p.upright = upright;
    drop(piece, x, z);
    if (!fits(piece)) {
        p = before;
        rebuild();
        return false;
    }
    return true;
}

void PieceLevel::unplace(int piece) {
    if (phase != Phase::Setup) return;
    cancel();
    Piece& p = pieces[static_cast<size_t>(piece)];
    p.placed = false;
    p.pos = p.home;
    p.yaw = 0;
    p.upright = false;
    rebuild();
}

int PieceLevel::piece_under(const Ray& ray) const {
    int best = -1;
    Real best_t = 1e30f;
    for (size_t i = 0; i < pieces.size(); ++i) {
        const Piece& p = pieces[i];
        const std::optional<RayHit> h = raycast(p.placed ? world.bodies[static_cast<size_t>(p.body)] : make_body(p, p.home), ray);
        if (h && h->t < best_t) {
            best_t = h->t;
            best = static_cast<int>(i);
        }
    }
    // Something else in front of it (a wall, another body) hides it.
    if (best >= 0) {
        const int own = pieces[static_cast<size_t>(best)].body;
        if (const auto first = raycast(world, ray, own); first && first->t < best_t - 1e-3f && first->body != 0) return -1;
    }
    return best;
}

void PieceLevel::press(const Ray& ray) {
    if (phase != Phase::Setup || held_ >= 0) return;
    const int i = piece_under(ray);
    if (i < 0) return;
    held_ = i;
    before_ = pieces[static_cast<size_t>(i)];
    last_spot_ = before_.pos;
    drag(ray);
}

void PieceLevel::drag(const Ray& ray) {
    if (held_ < 0) return;
    // The spot under the pointer: the first thing the ray meets other than the piece itself, or the ground.
    Piece& p = pieces[static_cast<size_t>(held_)];
    const std::optional<RayHit> hit = raycast(world, ray, p.body);
    if (hit) last_spot_ = hit->point;
    else if (const auto g = hit_plane_y(ray, kGroundY)) last_spot_ = *g;
    drop(held_, last_spot_.x, last_spot_.z);
    held_fits_ = fits(held_);
}

void PieceLevel::release(const Ray& ray) {
    if (held_ < 0) return;
    drag(ray);
    if (!held_fits_) {  // let go somewhere it cannot go: back where it came from
        cancel();
        return;
    }
    held_ = -1;
}

void PieceLevel::rotate(int direction) {
    if (held_ < 0) return;
    pieces[static_cast<size_t>(held_)].yaw += static_cast<Real>(direction) * rotation_step();
    drop(held_, last_spot_.x, last_spot_.z);
    held_fits_ = fits(held_);
}

void PieceLevel::tip() {
    if (held_ < 0) return;
    Piece& p = pieces[static_cast<size_t>(held_)];
    if (p.kind.radius > 0) return;
    p.upright = !p.upright;
    drop(held_, last_spot_.x, last_spot_.z);
    held_fits_ = fits(held_);
}

bool PieceLevel::cancel() {
    if (held_ < 0) return false;
    pieces[static_cast<size_t>(held_)] = before_;
    held_ = -1;
    rebuild();
    return true;
}

}  // namespace puzzle3d
