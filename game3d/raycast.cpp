#include "raycast.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace puzzle3d {

using phys3d::Body;
using phys3d::Real;
using phys3d::Shape;
using phys3d::Vec3;

namespace {

constexpr Real kInf = std::numeric_limits<Real>::infinity();

// Ray (in the shape's own frame) against a sphere at `c`. Returns the entry distance and normal.
bool sphere(Vec3 o, Vec3 d, Vec3 c, Real r, Real& t, Vec3& n) {
    const Vec3 m = o - c;
    const Real b = dot(m, d), cc = m.length_sq() - r * r;
    if (cc > 0 && b > 0) return false;
    const Real disc = b * b - cc;
    if (disc < 0) return false;
    t = std::max(Real(0), -b - std::sqrt(disc));
    n = (o + d * t - c).normalized();
    return true;
}

// Convex solid given as half-spaces dot(normal, p) <= offset: clip the ray's interval against each one.
template <class Planes>
bool planes(Vec3 o, Vec3 d, const Planes& for_each, Real& t, Vec3& n) {
    Real enter = 0, leave = kInf;
    Vec3 enter_normal;
    bool ok = true;
    for_each([&](Vec3 pn, Real off) {
        if (!ok) return;
        const Real denom = dot(pn, d), dist = off - dot(pn, o);
        if (std::fabs(denom) < 1e-9f) {
            if (dist < 0) ok = false;  // parallel and outside
            return;
        }
        const Real s = dist / denom;
        if (denom < 0) {
            if (s > enter) {
                enter = s;
                enter_normal = pn;
            }
        } else {
            leave = std::min(leave, s);
        }
        if (enter > leave) ok = false;
    });
    if (!ok) return false;
    t = enter;
    n = enter_normal;
    return true;
}

}  // namespace

std::optional<RayHit> raycast(const Body& body, const Ray& ray) {
    const Vec3 o = inv_rotate(body.q, ray.origin - body.pos), d = inv_rotate(body.q, ray.dir);
    Real t = 0;
    Vec3 n;
    bool hit = false;
    switch (body.shape.type) {
        case Shape::Type::Sphere: hit = sphere(o, d, {}, body.shape.radius, t, n); break;
        case Shape::Type::Box: {
            const Vec3 h = body.shape.half_extents;
            hit = planes(o, d, [&](auto&& f) {
                for (int axis = 0; axis < 3; ++axis)
                    for (const Real s : {Real(1), Real(-1)}) {
                        Vec3 pn;
                        pn[axis] = s;
                        f(pn, h[axis]);
                    }
            }, t, n);
            break;
        }
        case Shape::Type::Capsule: {
            // The two end spheres, and the cylinder between them (its side only: the ends are the spheres).
            const Real r = body.shape.radius, hl = body.shape.half_length;
            Real best = kInf;
            Vec3 best_n;
            for (const Real y : {hl, -hl}) {
                Real s;
                Vec3 sn;
                if (sphere(o, d, {0, y, 0}, r, s, sn) && s < best) {
                    best = s;
                    best_n = sn;
                }
            }
            const Real a = d.x * d.x + d.z * d.z, b = o.x * d.x + o.z * d.z, c = o.x * o.x + o.z * o.z - r * r;
            if (a > 1e-12f && b * b - a * c >= 0) {
                const Real s = std::max(Real(0), (-b - std::sqrt(b * b - a * c)) / a);
                const Vec3 p = o + d * s;
                if (std::fabs(p.y) <= hl && s < best && (c <= 0 || b < 0)) {
                    best = s;
                    best_n = Vec3{p.x, 0, p.z}.normalized();
                }
            }
            hit = best < kInf;
            t = best;
            n = best_n;
            break;
        }
        case Shape::Type::Hull:
            hit = planes(o, d, [&](auto&& f) {
                for (const phys3d::Hull::Face& face : body.shape.hull->faces) f(face.normal, face.offset);
            }, t, n);
            break;
    }
    if (!hit) return std::nullopt;
    RayHit out;
    out.t = t;
    out.point = ray.at(t);
    out.normal = rotate(body.q, n);
    return out;
}

std::optional<RayHit> raycast(const phys3d::World& world, const Ray& ray, int skip) {
    std::optional<RayHit> best;
    for (size_t i = 0; i < world.bodies.size(); ++i) {
        if (static_cast<int>(i) == skip) continue;
        std::optional<RayHit> h = raycast(world.bodies[i], ray);
        if (h && (!best || h->t < best->t)) {
            h->body = static_cast<int>(i);
            best = h;
        }
    }
    return best;
}

std::optional<Vec3> hit_plane_y(const Ray& ray, Real y) {
    if (std::fabs(ray.dir.y) < 1e-6f) return std::nullopt;
    const Real t = (y - ray.origin.y) / ray.dir.y;
    if (t < 0) return std::nullopt;
    return ray.at(t);
}

}  // namespace puzzle3d
