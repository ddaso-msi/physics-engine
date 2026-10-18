#include "painter.hpp"

#include <algorithm>
#include <cmath>

namespace puzzle3d::app {

using phys3d::Shape;

namespace {
const Vec3 kSun = Vec3{-0.45f, 1.0f, 0.35f}.normalized();
}

Painter::Painter(SDL_Renderer* ren, const demo3d::Camera& cam) : ren_(ren), cam_(cam), eye_(cam.eye()) {}

bool Painter::project(Vec3 p, SDL_FPoint& out, Real* depth) const {
    const Vec3 v = cam_.to_view(p);
    if (v.z < cam_.near_plane) return false;
    const demo3d::ScreenPoint s = cam_.to_screen(v);
    out = {s.x, s.y};
    if (depth) *depth = v.z;
    return true;
}

SDL_FColor Painter::lit(SDL_Color c, Uint8 alpha, Vec3 n) const {
    const float light = 0.5f + 0.5f * std::max(0.0f, dot(n, kSun)) + 0.08f * n.y;
    auto ch = [&](Uint8 v) { return std::min(1.0f, static_cast<float>(v) / 255.0f * light); };
    return {ch(c.r), ch(c.g), ch(c.b), static_cast<float>(alpha) / 255.0f};
}

void Painter::tri(Vec3 a, Vec3 b, Vec3 c, Vec3 na, Vec3 nb, Vec3 nc, Shade s) {
    Tri t;
    Real za, zb, zc;
    SDL_FPoint pa, pb, pc;
    if (!project(a, pa, &za) || !project(b, pb, &zb) || !project(c, pc, &zc)) return;
    t.v[0] = {pa, lit(s.base, s.alpha, na), {0, 0}};
    t.v[1] = {pb, lit(s.base, s.alpha, nb), {0, 0}};
    t.v[2] = {pc, lit(s.base, s.alpha, nc), {0, 0}};
    t.depth = (za + zb + zc) / 3.0f;
    tris_.push_back(t);
}

void Painter::quad(Vec3 a, Vec3 b, Vec3 c, Vec3 d, Vec3 n, Shade s, int tu, int tv) {
    // a b c d go round the face; cut it into tu x tv tiles so a big face sorts by its parts.
    if (dot(n, eye_ - (a + c) * 0.5f) <= 0 && tu == 1 && tv == 1) return;  // turned away
    for (int i = 0; i < tu; ++i)
        for (int k = 0; k < tv; ++k) {
            auto at = [&](Real u, Real v) { return a + (b - a) * u + (d - a) * v; };
            const Real u0 = static_cast<Real>(i) / static_cast<Real>(tu), u1 = static_cast<Real>(i + 1) / static_cast<Real>(tu);
            const Real v0 = static_cast<Real>(k) / static_cast<Real>(tv), v1 = static_cast<Real>(k + 1) / static_cast<Real>(tv);
            const Vec3 p00 = at(u0, v0), p10 = at(u1, v0), p11 = at(u1, v1), p01 = at(u0, v1);
            if (dot(n, eye_ - (p00 + p11) * 0.5f) <= 0) continue;
            tri(p00, p10, p11, n, n, n, s);
            tri(p00, p11, p01, n, n, n, s);
        }
}

void Painter::sky() {
    SDL_Vertex v[4];
    const SDL_FColor top{color::kSky.r / 255.0f, color::kSky.g / 255.0f, color::kSky.b / 255.0f, 1};
    const SDL_FColor low{color::kSkyLow.r / 255.0f, color::kSkyLow.g / 255.0f, color::kSkyLow.b / 255.0f, 1};
    const float w = static_cast<float>(kWidth), h = static_cast<float>(kHeight);
    v[0] = {{0, 0}, top, {0, 0}};
    v[1] = {{w, 0}, top, {0, 0}};
    v[2] = {{w, h}, low, {0, 0}};
    v[3] = {{0, h}, low, {0, 0}};
    const int idx[6] = {0, 1, 2, 0, 2, 3};
    SDL_RenderGeometry(ren_, nullptr, v, 4, idx, 6);
}

void Painter::ground(const Body& g) {
    // A checkerboard of metre tiles, fading toward the edge of the slab. Drawn straight away: nothing is under it.
    const Real y = g.pos.y + g.shape.half_extents.y;
    const int reach = static_cast<int>(g.shape.half_extents.x);
    std::vector<SDL_Vertex> verts;
    std::vector<int> idx;
    for (int x = -reach; x < reach; ++x)
        for (int z = -reach; z < reach; ++z) {
            const Vec3 c[4] = {{static_cast<Real>(x), y, static_cast<Real>(z)}, {static_cast<Real>(x + 1), y, static_cast<Real>(z)},
                               {static_cast<Real>(x + 1), y, static_cast<Real>(z + 1)}, {static_cast<Real>(x), y, static_cast<Real>(z + 1)}};
            SDL_FPoint p[4];
            bool ok = true;
            for (int k = 0; k < 4 && ok; ++k) ok = project(c[k], p[k]);
            if (!ok) continue;
            const float fade = std::max(0.0f, 1.0f - std::sqrt(static_cast<float>(x * x + z * z)) / static_cast<float>(reach));
            const bool odd = ((x + z) & 1) != 0;
            const float base = (odd ? 0.27f : 0.30f) * (0.55f + 0.45f * fade);
            const SDL_FColor col{base * 0.92f, base, base * 1.12f, 1};
            const int first = static_cast<int>(verts.size());
            for (int k = 0; k < 4; ++k) verts.push_back({p[k], col, {0, 0}});
            for (int k : {0, 1, 2, 0, 2, 3}) idx.push_back(first + k);
        }
    SDL_RenderGeometry(ren_, nullptr, verts.data(), static_cast<int>(verts.size()), idx.data(), static_cast<int>(idx.size()));
}

void Painter::shadow(const Body& b) {
    const Real r = b.bounding_radius() * 0.85f, height = b.pos.y - kGroundY;
    if (height > 15 || height < -0.5f) return;
    const float strength = std::clamp(0.45f * (1.0f - static_cast<float>(height) / 15.0f), 0.0f, 0.45f);
    constexpr int kSegments = 20;
    SDL_Vertex v[kSegments + 1];
    int idx[kSegments * 3];
    SDL_FPoint centre;
    if (!project({b.pos.x, kGroundY + 0.003f, b.pos.z}, centre)) return;
    v[0] = {centre, {0, 0, 0, strength}, {0, 0}};
    for (int i = 0; i < kSegments; ++i) {
        const Real a = 6.2831853f * static_cast<Real>(i) / kSegments;
        SDL_FPoint p;
        if (!project({b.pos.x + r * std::cos(a), kGroundY + 0.003f, b.pos.z + r * std::sin(a)}, p)) return;
        v[i + 1] = {p, {0, 0, 0, 0}, {0, 0}};
        idx[i * 3] = 0;
        idx[i * 3 + 1] = i + 1;
        idx[i * 3 + 2] = (i + 1) % kSegments + 1;
    }
    SDL_RenderGeometry(ren_, nullptr, v, kSegments + 1, idx, kSegments * 3);
}

void Painter::round_body(const Body& b, Shade s) {
    // Rings of latitude. A capsule is a sphere cut at its equator, the two halves pushed apart by the half
    // length, with a cylinder between them.
    constexpr int kAround = 18, kRings = 10;
    const Real r = b.shape.radius, hl = b.shape.type == Shape::Type::Capsule ? b.shape.half_length : Real(0);
    const phys3d::Transform t = b.transform();
    auto ring_point = [&](int ring, int k, Real shift, Vec3& normal) {
        const Real phi = 3.14159265f * (static_cast<Real>(ring) / kRings - 0.5f);  // -pi/2 .. pi/2
        const Real th = 6.2831853f * static_cast<Real>(k) / kAround;
        const Vec3 n{std::cos(phi) * std::cos(th), std::sin(phi), std::cos(phi) * std::sin(th)};
        normal = rotate(b.q, n);
        return apply(t, n * r + Vec3{0, shift, 0});
    };
    for (const int half : {0, 1}) {
        const Real shift = half == 0 ? -hl : hl;
        for (int ring = half * kRings / 2; ring < (half + 1) * kRings / 2; ++ring)
            for (int k = 0; k < kAround; ++k) {
                Vec3 n00, n10, n11, n01;
                const Vec3 p00 = ring_point(ring, k, shift, n00), p10 = ring_point(ring, k + 1, shift, n10);
                const Vec3 p11 = ring_point(ring + 1, k + 1, shift, n11), p01 = ring_point(ring + 1, k, shift, n01);
                const Vec3 mid_n = (n00 + n11).normalized();
                if (dot(mid_n, eye_ - (p00 + p11) * 0.5f) <= 0) continue;
                tri(p00, p11, p10, n00, n11, n10, s);
                tri(p00, p01, p11, n00, n01, n11, s);
            }
    }
    if (hl <= 0) return;
    for (int k = 0; k < kAround; ++k) {
        Vec3 w0, w1, unused;
        const Vec3 lo0 = ring_point(kRings / 2, k, -hl, w0), lo1 = ring_point(kRings / 2, k + 1, -hl, w1);
        const Vec3 hi0 = ring_point(kRings / 2, k, hl, unused), hi1 = ring_point(kRings / 2, k + 1, hl, unused);
        if (dot((w0 + w1).normalized(), eye_ - (lo0 + hi1) * 0.5f) <= 0) continue;
        tri(lo0, hi1, lo1, w0, w1, w1, s);
        tri(lo0, hi0, hi1, w0, w0, w1, s);
    }
}

void Painter::body(const Body& b, Shade s) {
    switch (b.shape.type) {
        case Shape::Type::Box: {
            const Vec3 h = b.shape.half_extents;
            const phys3d::Transform t = b.transform();
            for (int axis = 0; axis < 3; ++axis)
                for (const Real sign : {Real(1), Real(-1)}) {
                    const int u = (axis + 1) % 3, v = (axis + 2) % 3;
                    Vec3 n, du, dv, centre;
                    n[axis] = sign;
                    centre[axis] = sign * h[axis];
                    du[u] = h[u];
                    dv[v] = h[v] * sign;  // keeps the corners going round the same way seen from outside
                    const Vec3 corners[4] = {apply(t, centre - du - dv), apply(t, centre + du - dv), apply(t, centre + du + dv), apply(t, centre - du + dv)};
                    const int tu = std::max(1, static_cast<int>(h[u] * 2 / 1.0f)), tv = std::max(1, static_cast<int>(h[v] * 2 / 1.0f));
                    quad(corners[0], corners[1], corners[2], corners[3], rotate(b.q, n), s, tu, tv);
                }
            break;
        }
        case Shape::Type::Sphere:
        case Shape::Type::Capsule: round_body(b, s); break;
        case Shape::Type::Hull: {
            const phys3d::Hull& hull = *b.shape.hull;
            for (const phys3d::Hull::Face& f : hull.faces) {
                const Vec3 n = rotate(b.q, f.normal);
                const Vec3 p0 = apply(b.transform(), hull.point(f, 0));
                if (dot(n, eye_ - p0) <= 0) continue;
                for (int i = 1; i + 1 < f.count; ++i)
                    tri(p0, apply(b.transform(), hull.point(f, i)), apply(b.transform(), hull.point(f, i + 1)), n, n, n, s);
            }
            break;
        }
    }
}

void Painter::zone_floor(Vec3 lo, Vec3 hi, SDL_Color c) {
    const Real y = lo.y + 0.01f;
    quad({lo.x, y, lo.z}, {lo.x, y, hi.z}, {hi.x, y, hi.z}, {hi.x, y, lo.z}, {0, 1, 0}, {c, 70}, 3, 3);
}

void Painter::flush() {
    std::sort(tris_.begin(), tris_.end(), [](const Tri& a, const Tri& b) { return a.depth > b.depth; });
    std::vector<SDL_Vertex> verts;
    verts.reserve(tris_.size() * 3);
    for (const Tri& t : tris_)
        for (const SDL_Vertex& v : t.v) verts.push_back(v);
    if (!verts.empty()) SDL_RenderGeometry(ren_, nullptr, verts.data(), static_cast<int>(verts.size()), nullptr, 0);
    tris_.clear();
}

void Painter::line(Vec3 a, Vec3 b, SDL_Color c, bool dashed, float width) {
    const Real length = distance(a, b);
    const int pieces = dashed ? std::max(1, static_cast<int>(length / 0.18f)) : 1;
    set_color(ren_, c);
    for (int i = 0; i < pieces; i += dashed ? 2 : 1) {
        const Vec3 p = a + (b - a) * (static_cast<Real>(i) / static_cast<Real>(pieces));
        const Vec3 q = a + (b - a) * (static_cast<Real>(i + 1) / static_cast<Real>(pieces));
        demo3d::ScreenPoint sa, sb;
        if (!cam_.project_segment(p, q, sa, sb)) continue;
        const float dx = sb.x - sa.x, dy = sb.y - sa.y, len = std::sqrt(dx * dx + dy * dy);
        if (len < 0.01f) continue;
        for (float o = -(width - 1) * 0.5f; o <= (width - 1) * 0.5f + 0.01f; o += 1.0f)
            SDL_RenderLine(ren_, sa.x - dy / len * o, sa.y + dx / len * o, sb.x - dy / len * o, sb.y + dx / len * o);
    }
}

void Painter::marker(Vec3 p, Real radius, SDL_Color c) {
    SDL_FPoint s;
    Real z;
    if (!project(p, s, &z)) return;
    fill_disc(ren_, s, std::max(2.0f, static_cast<float>(radius * cam_.focal() / z)), c);
}

void Painter::label(Vec3 p, const char* str, SDL_Color c) {
    SDL_FPoint s;
    if (!project(p, s)) return;
    const float w = text_width(str, 2.0f);
    fill_rect(ren_, {s.x - w * 0.5f - 4, s.y - 10, w + 8, 20}, {18, 21, 30, 170});
    text_centered(ren_, s.x, s.y - 8, 2.0f, c, str);
}

}  // namespace puzzle3d::app
