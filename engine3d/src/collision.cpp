#include <phys3d/collision.hpp>

#include <algorithm>
#include <limits>

namespace phys3d {

namespace {

constexpr Real kRealMax = std::numeric_limits<Real>::max();

// Preference between separating-axis candidates. A face axis is kept unless another candidate is
// CLEARLY better (5% and 1 cm less overlap). Face contacts give a whole patch of points and a normal that
// does not swing about; without the bias, two nearly tied axes would trade places from frame to frame
// and the contact would flicker between a 4-point patch and a single edge point.
constexpr Real kRelativeTolerance = static_cast<Real>(0.95);
constexpr Real kAbsoluteTolerance = static_cast<Real>(0.01);

// The side planes of a reference face are pushed out by this much before clipping. When two equal boxes
// are stacked, the corners of one face lie exactly ON the side planes of the other, and rounding would
// put each corner just inside one frame (kept as "that corner") and just outside the next (replaced by "a
// point cut by plane k"), a different feature with a different id each time. With the margin a corner
// that is within 2 mm of the edge simply stays a corner.
constexpr Real kClipMargin = static_cast<Real>(0.002);

Real sign_of(Real x) { return x < 0 ? Real(-1) : Real(1); }

// ---------------------------------------------------------------------------------------------
// sphere vs sphere
// ---------------------------------------------------------------------------------------------
bool collide_spheres(const Body& a, const Body& b, Manifold& out) {
    const Real ra = a.shape.radius, rb = b.shape.radius, r = ra + rb;
    const Vec3 d = b.pos - a.pos;
    const Real dist_sq = d.length_sq();
    if (dist_sq > r * r) return false;
    const Real dist = std::sqrt(dist_sq);
    const Vec3 n = dist > static_cast<Real>(1e-9) ? d / dist : Vec3{0, 1, 0};  // concentric: any way out will do
    out.normal = n;
    out.depth = r - dist;
    out.count = 1;
    out.points[0] = {a.pos + n * (ra - out.depth * static_cast<Real>(0.5)), out.depth, 0};
    return true;
}

// ---------------------------------------------------------------------------------------------
// box vs sphere (normal points from the box to the sphere)
// ---------------------------------------------------------------------------------------------
bool collide_box_sphere(const Body& box, const Body& sphere, Manifold& out) {
    const Vec3 h = box.shape.half_extents;
    const Real r = sphere.shape.radius;
    // In the box's own frame it is an axis-aligned block, and the nearest point of the block to the
    // sphere's centre is just the centre with each coordinate clamped to the block.
    const Vec3 c = inv_rotate(box.q, sphere.pos - box.pos);
    const Vec3 closest{std::clamp(c.x, -h.x, h.x), std::clamp(c.y, -h.y, h.y), std::clamp(c.z, -h.z, h.z)};
    const Vec3 d = c - closest;
    const Real dist_sq = d.length_sq();
    if (dist_sq > r * r) return false;

    Vec3 n_local;
    Real s;  // signed distance from the box surface to the sphere centre, along n_local
    std::uint32_t id = 0;
    if (dist_sq > static_cast<Real>(1e-12)) {
        // Centre outside: the nearest feature is a face, an edge or a corner depending on how many
        // coordinates had to be clamped (1, 2 or 3). The id records which, and on which side.
        s = std::sqrt(dist_sq);
        n_local = d / s;
        for (int i = 0; i < 3; ++i) id |= static_cast<std::uint32_t>(c[i] > h[i] ? 2 : c[i] < -h[i] ? 1 : 0) << (2 * i);
    } else {
        // Centre inside: leave through the nearest face.
        int axis = 0;
        Real nearest = kRealMax;
        for (int i = 0; i < 3; ++i) {
            const Real to_face = h[i] - std::fabs(c[i]);
            if (to_face < nearest) {
                nearest = to_face;
                axis = i;
            }
        }
        n_local = {};
        n_local[axis] = sign_of(c[axis]);
        s = -nearest;
        id = 0x100u | static_cast<std::uint32_t>(axis << 1) | (c[axis] < 0 ? 1u : 0u);
    }

    const Vec3 n = rotate(box.q, n_local);
    out.normal = n;
    out.depth = r - s;
    out.count = 1;
    // The sphere's deepest point is centre - n r, the box surface is centre - n s; take the midpoint.
    out.points[0] = {sphere.pos - n * ((r + s) * static_cast<Real>(0.5)), out.depth, id};
    return true;
}

// ---------------------------------------------------------------------------------------------
// box vs box
// ---------------------------------------------------------------------------------------------
struct WorldBox {
    Vec3 center;
    Vec3 axis[3];  // unit vectors: the box's own x, y, z in world space
    Vec3 half;
};

WorldBox to_world(const Body& b) {
    const Mat3 r = to_mat3(b.q);
    return {b.pos, {r.cx, r.cy, r.cz}, b.shape.half_extents};
}

// Half the length of the box's shadow on unit axis `l`.
Real projected_radius(const WorldBox& box, Vec3 l) {
    return box.half.x * std::fabs(dot(box.axis[0], l)) + box.half.y * std::fabs(dot(box.axis[1], l)) +
           box.half.z * std::fabs(dot(box.axis[2], l));
}

// A point of a polygon being clipped, with a tag saying where it came from (an original corner of the
// incident face, or which plane cut which edge), so the same contact point keeps the same id next frame.
struct ClipVertex {
    Vec3 p;
    std::uint32_t tag = 0;
};

// Sutherland-Hodgman against one plane: keep the part of the polygon where dot(n, p) <= c.
int clip_polygon(ClipVertex* out, const ClipVertex* in, int count, Vec3 n, Real c, std::uint32_t plane) {
    int kept = 0;
    for (int i = 0; i < count; ++i) {
        const ClipVertex& cur = in[i];
        const ClipVertex& next = in[(i + 1) % count];
        const Real d_cur = dot(n, cur.p) - c, d_next = dot(n, next.p) - c;
        if (d_cur <= 0) out[kept++] = cur;
        if ((d_cur < 0 && d_next > 0) || (d_cur > 0 && d_next < 0)) {  // the edge crosses the plane
            const Real t = d_cur / (d_cur - d_next);
            out[kept++] = {cur.p + (next.p - cur.p) * t, 1000u + (cur.tag * 37u + next.tag) * 37u + plane};
        }
    }
    return kept;
}

// Up to 8 candidate points come out of clipping; keep the 4 that preserve the most of the contact patch:
// the deepest, the one farthest from it, the one that makes the biggest triangle with those two, and
// the one farthest out on the other side of that first edge.
int reduce_to_four(ContactPoint* pts, int count, Vec3 u, Vec3 v) {
    if (count <= Manifold::kMaxPoints) return count;
    auto flat = [&](int i) { return phys::Vec2{dot(pts[i].point, u), dot(pts[i].point, v)}; };
    bool used[16] = {};
    int chosen[4];

    chosen[0] = 0;
    for (int i = 1; i < count; ++i)
        if (pts[i].depth > pts[chosen[0]].depth) chosen[0] = i;
    used[chosen[0]] = true;

    chosen[1] = -1;
    Real best = -1;
    for (int i = 0; i < count; ++i) {
        if (used[i]) continue;
        const Real d = (flat(i) - flat(chosen[0])).length_sq();
        if (d > best) {
            best = d;
            chosen[1] = i;
        }
    }
    used[chosen[1]] = true;

    const phys::Vec2 a = flat(chosen[0]), edge = flat(chosen[1]) - a;
    auto side = [&](int i) { return phys::cross(edge, flat(i) - a); };  // twice the signed triangle area
    chosen[2] = -1;
    best = -1;
    for (int i = 0; i < count; ++i) {
        if (used[i]) continue;
        if (std::fabs(side(i)) > best) {
            best = std::fabs(side(i));
            chosen[2] = i;
        }
    }
    used[chosen[2]] = true;
    const Real third_side = side(chosen[2]);

    chosen[3] = -1;
    best = 0;
    for (int i = 0; i < count; ++i) {  // farthest on the opposite side of the first edge
        if (used[i]) continue;
        const Real opposite = third_side >= 0 ? -side(i) : side(i);
        if (opposite > best) {
            best = opposite;
            chosen[3] = i;
        }
    }
    if (chosen[3] < 0) {  // everything left is on the same side: take the one farthest from the triangle's corners
        best = -1;
        for (int i = 0; i < count; ++i) {
            if (used[i]) continue;
            Real nearest = kRealMax;
            for (int k = 0; k < 3; ++k) nearest = std::min(nearest, (flat(i) - flat(chosen[k])).length_sq());
            if (nearest > best) {
                best = nearest;
                chosen[3] = i;
            }
        }
    }

    ContactPoint kept[4];
    for (int k = 0; k < 4; ++k) kept[k] = pts[chosen[k]];
    for (int k = 0; k < 4; ++k) pts[k] = kept[k];
    return 4;
}

// Face contact: `ref`'s face number `ref_axis` (on the side facing `inc`) is the reference face. Clip the
// most opposed face of `inc` against the four side planes of the reference face; what is left and lies
// behind the reference face is the contact patch. `n` is the reference face's outward normal.
void face_manifold(const WorldBox& ref, const WorldBox& inc, int ref_axis, Vec3 n, Real overlap, bool ref_is_a,
                   Manifold& out) {
    // The incident face: the face of `inc` whose normal points most directly against n.
    int inc_axis = 0;
    Real most = -1;
    for (int k = 0; k < 3; ++k) {
        const Real d = std::fabs(dot(inc.axis[k], n));
        if (d > most) {
            most = d;
            inc_axis = k;
        }
    }
    const Real inc_sign = dot(inc.axis[inc_axis], n) > 0 ? Real(-1) : Real(1);
    const int iu = (inc_axis + 1) % 3, iv = (inc_axis + 2) % 3;
    const Vec3 inc_center = inc.center + inc.axis[inc_axis] * (inc_sign * inc.half[inc_axis]);
    const Vec3 eu = inc.axis[iu] * inc.half[iu], ev = inc.axis[iv] * inc.half[iv];

    ClipVertex poly_a[16], poly_b[16];
    poly_a[0] = {inc_center + eu + ev, 0};
    poly_a[1] = {inc_center - eu + ev, 1};
    poly_a[2] = {inc_center - eu - ev, 2};
    poly_a[3] = {inc_center + eu - ev, 3};
    int count = 4;

    // The reference face is a rectangle centred on ref_center spanned by its two in-plane axes.
    const int ru = (ref_axis + 1) % 3, rv = (ref_axis + 2) % 3;
    const Vec3 ref_center = ref.center + n * ref.half[ref_axis];
    const Vec3 side_normal[4] = {ref.axis[ru], -ref.axis[ru], ref.axis[rv], -ref.axis[rv]};
    const Real side_half[4] = {ref.half[ru], ref.half[ru], ref.half[rv], ref.half[rv]};
    ClipVertex *src = poly_a, *dst = poly_b;
    for (std::uint32_t k = 0; k < 4 && count > 0; ++k) {
        count = clip_polygon(dst, src, count, side_normal[k], dot(side_normal[k], ref_center) + side_half[k] + kClipMargin, k);
        std::swap(src, dst);
    }

    const std::uint32_t pair_code = (ref_is_a ? 0u : 1u << 30) | (static_cast<std::uint32_t>(ref_axis) << 26) |
                                    (dot(ref.axis[ref_axis], n) > 0 ? 0u : 1u << 25) |
                                    (static_cast<std::uint32_t>(inc_axis) << 22) | (inc_sign > 0 ? 0u : 1u << 21);

    ContactPoint candidates[16];
    int found = 0;
    for (int i = 0; i < count; ++i) {
        const Real dist = dot(n, src[i].p - ref_center);  // <= 0 means behind the reference face
        if (dist <= 0)
            candidates[found++] = {src[i].p - n * (dist * static_cast<Real>(0.5)), -dist,
                                   pair_code | (src[i].tag & 0x1FFFFFu)};
    }
    if (found == 0) {
        // SAT proved an overlap, but the incident face's patch over the reference face is empty (the boxes
        // only clip corners). Fall back to the deepest corner of the incident face.
        int deepest = 0;
        for (int i = 1; i < 4; ++i)
            if (dot(n, poly_a[i].p) < dot(n, poly_a[deepest].p)) deepest = i;
        const Real dist = std::min(dot(n, poly_a[deepest].p - ref_center), Real(0));
        candidates[found++] = {poly_a[deepest].p - n * (dist * static_cast<Real>(0.5)), overlap,
                               pair_code | 0x100000u | static_cast<std::uint32_t>(deepest)};
    }
    found = reduce_to_four(candidates, found, ref.axis[ru], ref.axis[rv]);

    out.normal = ref_is_a ? n : -n;  // always A -> B
    out.depth = overlap;
    out.count = found;
    for (int i = 0; i < found; ++i) out.points[i] = candidates[i];
}

// Edge contact: an edge of A (parallel to its axis ia) against an edge of B (parallel to its axis ib).
// The contact is the midpoint of the shortest segment between the two edge lines.
void edge_manifold(const WorldBox& a, const WorldBox& b, int ia, int ib, Vec3 n, Real overlap, Manifold& out) {
    // The edge of A farthest along n, and the edge of B farthest against it.
    Vec3 pa = a.center, pb = b.center;
    std::uint32_t id = 0x80000000u | (static_cast<std::uint32_t>(ia) << 8) | (static_cast<std::uint32_t>(ib) << 4);
    for (int k = 0; k < 3; ++k) {
        if (k != ia) {
            const Real s = sign_of(dot(a.axis[k], n));
            pa += a.axis[k] * (s * a.half[k]);
            if (s < 0) id |= 1u << (12 + k);
        }
        if (k != ib) {
            const Real s = sign_of(dot(b.axis[k], n));
            pb -= b.axis[k] * (s * b.half[k]);
            if (s < 0) id |= 1u << (16 + k);
        }
    }
    // Closest points of the lines pa + s da and pb + t db (da, db unit, not parallel).
    const Vec3 da = a.axis[ia], db = b.axis[ib], r = pa - pb;
    const Real bb = dot(da, db), c = dot(da, r), f = dot(db, r);
    const Real denom = 1 - bb * bb;
    Real s = denom > static_cast<Real>(1e-9) ? (bb * f - c) / denom : 0;
    Real t = denom > static_cast<Real>(1e-9) ? (f - bb * c) / denom : 0;
    s = std::clamp(s, -a.half[ia], a.half[ia]);
    t = std::clamp(t, -b.half[ib], b.half[ib]);

    out.normal = n;
    out.depth = overlap;
    out.count = 1;
    out.points[0] = {(pa + da * s + pb + db * t) * static_cast<Real>(0.5), overlap, id};
}

bool collide_boxes(const Body& body_a, const Body& body_b, Manifold& out) {
    const WorldBox a = to_world(body_a), b = to_world(body_b);
    const Vec3 t = b.center - a.center;

    // Separating axis test. Two convex boxes are apart exactly when their shadows on some axis do not
    // overlap, and only 15 axes need checking: the 3 face normals of each box, and the 9 directions
    // perpendicular to one edge of each (the face normals alone miss boxes that pass edge to edge).
    // `overlap` on an axis is how far the shadows interpenetrate; the smallest one is the way out.
    Real face_a = kRealMax, face_b = kRealMax, edge = kRealMax;
    int axis_a = 0, axis_b = 0, edge_ia = 0, edge_ib = 0;
    Vec3 edge_axis;

    for (int i = 0; i < 3; ++i) {
        const Real overlap = a.half[i] + projected_radius(b, a.axis[i]) - std::fabs(dot(t, a.axis[i]));
        if (overlap < 0) return false;
        if (overlap < face_a) {
            face_a = overlap;
            axis_a = i;
        }
    }
    for (int j = 0; j < 3; ++j) {
        const Real overlap = projected_radius(a, b.axis[j]) + b.half[j] - std::fabs(dot(t, b.axis[j]));
        if (overlap < 0) return false;
        if (overlap < face_b) {
            face_b = overlap;
            axis_b = j;
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Vec3 l = cross(a.axis[i], b.axis[j]);
            const Real len_sq = l.length_sq();
            if (len_sq < static_cast<Real>(1e-6)) continue;  // parallel edges: a face axis already covers this
            l = l / std::sqrt(len_sq);
            const Real overlap = projected_radius(a, l) + projected_radius(b, l) - std::fabs(dot(t, l));
            if (overlap < 0) return false;
            if (overlap < edge) {
                edge = overlap;
                edge_ia = i;
                edge_ib = j;
                edge_axis = l;
            }
        }
    }

    // Choose the contact type: A's face unless B's face, and then an edge pair, is clearly shallower.
    enum class Kind { FaceA, FaceB, Edge } kind = Kind::FaceA;
    Real best = face_a;
    if (face_b < kRelativeTolerance * best - kAbsoluteTolerance) {
        kind = Kind::FaceB;
        best = face_b;
    }
    if (edge < kRelativeTolerance * best - kAbsoluteTolerance) {
        kind = Kind::Edge;
        best = edge;
    }

    Manifold m;
    if (kind == Kind::FaceA) {
        const Vec3 n = a.axis[axis_a] * sign_of(dot(a.axis[axis_a], t));  // A's face that looks at B
        face_manifold(a, b, axis_a, n, best, true, m);
    } else if (kind == Kind::FaceB) {
        const Vec3 n = b.axis[axis_b] * -sign_of(dot(b.axis[axis_b], t));  // B's face that looks at A
        face_manifold(b, a, axis_b, n, best, false, m);
    } else {
        edge_manifold(a, b, edge_ia, edge_ib, edge_axis * sign_of(dot(edge_axis, t)), best, m);
    }
    out = m;
    return true;
}

}  // namespace

bool collide(const Body& a, const Body& b, Manifold& out) {
    using Type = Shape::Type;
    Manifold m;
    bool hit;
    if (a.shape.type == Type::Sphere && b.shape.type == Type::Sphere) {
        hit = collide_spheres(a, b, m);
    } else if (a.shape.type == Type::Box && b.shape.type == Type::Sphere) {
        hit = collide_box_sphere(a, b, m);
    } else if (a.shape.type == Type::Sphere && b.shape.type == Type::Box) {
        hit = collide_box_sphere(b, a, m);
        m.normal = -m.normal;  // computed box -> sphere; we promise A -> B
    } else {
        hit = collide_boxes(a, b, m);
    }
    if (hit) out = m;
    return hit;
}

bool contains(const Body& body, Vec3 world_point) {
    const Vec3 local = inv_rotate(body.q, world_point - body.pos);
    if (body.shape.type == Shape::Type::Sphere) return local.length_sq() <= body.shape.radius * body.shape.radius;
    const Vec3 h = body.shape.half_extents;
    return std::fabs(local.x) <= h.x && std::fabs(local.y) <= h.y && std::fabs(local.z) <= h.z;
}

}  // namespace phys3d
