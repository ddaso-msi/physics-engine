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
    // Every "pick the best" loop below compares with >, which is false for NaN. If a body's position has
    // gone NaN (a blown-up simulation) no candidate would win, so fall back to any point not yet taken
    // rather than indexing with -1.
    auto first_unused = [&] {
        for (int i = 0; i < count; ++i)
            if (!used[i]) return i;
        return 0;
    };

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
    if (chosen[1] < 0) chosen[1] = first_unused();
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
    if (chosen[2] < 0) chosen[2] = first_unused();
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

    if (chosen[3] < 0) chosen[3] = first_unused();

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

// ---------------------------------------------------------------------------------------------
// any two convex shapes, by GJK and EPA
// ---------------------------------------------------------------------------------------------

// The flat piece of a shape's core that faces a given direction: a point (a sphere's centre, a hull's
// farthest vertex), a segment (a capsule's axis) or a rectangle (the face of a box that faces that way most
// squarely).
struct Feature {
    Vec3 p[4];
    int count = 0;  // 1, 2 or 4 (corners in order round the face)
    // For a face:
    Vec3 center, normal, u, v;  // normal is outward; u and v lie in the face
    Real hu = 0, hv = 0;        // half the face's size along u and v
    Real alignment = 0;         // how squarely it faces the direction asked for: 1 is exactly
    std::uint32_t code = 7;     // which face of the box (axis and side), for contact ids; 7 = not a face
};

Feature facing_feature(const Convex& c, Vec3 d) {
    Feature f;
    if (c.core == Convex::Core::Segment) {
        const Vec3 h = rotate(c.q, {0, c.half.y, 0});
        f.p[0] = c.pos - h;
        f.p[1] = c.pos + h;
        f.count = 2;
    } else if (c.core == Convex::Core::Box) {
        const Mat3 r = to_mat3(c.q);
        const Vec3 axis[3] = {r.cx, r.cy, r.cz};
        int k = 0;
        for (int i = 1; i < 3; ++i)
            if (std::fabs(dot(axis[i], d)) > std::fabs(dot(axis[k], d))) k = i;
        const bool negative = dot(axis[k], d) < 0;
        f.normal = negative ? -axis[k] : axis[k];
        f.alignment = std::fabs(dot(axis[k], d));
        f.u = axis[(k + 1) % 3];
        f.v = axis[(k + 2) % 3];
        f.hu = c.half[(k + 1) % 3];
        f.hv = c.half[(k + 2) % 3];
        f.center = c.pos + f.normal * c.half[k];
        const Vec3 eu = f.u * f.hu, ev = f.v * f.hv;
        f.p[0] = f.center + eu + ev;
        f.p[1] = f.center - eu + ev;
        f.p[2] = f.center - eu - ev;
        f.p[3] = f.center + eu - ev;
        f.count = 4;
        f.code = static_cast<std::uint32_t>(k) * 2 + (negative ? 1u : 0u);
    } else {
        // A point core is its own feature. A hull is only a cloud of points here, with no record of its
        // faces, so it offers its farthest vertex and gets a one-point contact.
        f.p[0] = c.core == Convex::Core::Point ? c.pos : c.core_support(d);
        f.count = 1;
    }
    return f;
}

// A face counts as "the" contact feature when the contact normal is within about 0.8 degrees of its own.
// When a face really is what is touching, EPA's normal IS that face's normal to rounding, so this need not
// be generous; and it must not be, because the manifold takes the face's normal in place of EPA's.
constexpr Real kFaceAligned = static_cast<Real>(0.9999);

// Keeps the part of segment (a, b) where dot(n, p) <= c. False if nothing is left.
bool clip_segment(ClipVertex& a, ClipVertex& b, Vec3 n, Real c, std::uint32_t plane) {
    const Real da = dot(n, a.p) - c, db = dot(n, b.p) - c;
    if (da > 0 && db > 0) return false;
    if (da > 0) a = {a.p + (b.p - a.p) * (da / (da - db)), 500u + a.tag * 8u + plane};
    if (db > 0) b = {b.p + (a.p - b.p) * (db / (db - da)), 500u + b.tag * 8u + plane};
    return true;
}

// The contact patch when one shape presents a face: clip the other shape's feature to the face's outline
// and keep what has sunk below it. This is the box code's face_manifold with the incident "face" allowed to
// be a segment or a single point, and with both shapes' rounding radii.
int face_patch(const Feature& ref, Real ref_radius, const Feature& inc, Real inc_radius, bool ref_is_a, ContactPoint* out) {
    ClipVertex buffer_a[16], buffer_b[16];
    ClipVertex *src = buffer_a, *dst = buffer_b;
    int count = inc.count;
    for (int i = 0; i < count; ++i) src[i] = {inc.p[i], static_cast<std::uint32_t>(i)};

    const Vec3 side_normal[4] = {ref.u, -ref.u, ref.v, -ref.v};
    const Real side_half[4] = {ref.hu, ref.hu, ref.hv, ref.hv};
    for (std::uint32_t k = 0; k < 4 && count > 0; ++k) {
        const Real limit = dot(side_normal[k], ref.center) + side_half[k] + kClipMargin;
        if (count == 2) {
            if (!clip_segment(src[0], src[1], side_normal[k], limit, k)) count = 0;
        } else if (count == 1) {
            if (dot(side_normal[k], src[0].p) > limit) count = 0;
        } else {
            count = clip_polygon(dst, src, count, side_normal[k], limit, k);
            std::swap(src, dst);
        }
    }

    const std::uint32_t pair_code = (ref_is_a ? 0u : 1u << 30) | (ref.code << 25) | (inc.code << 21);
    ContactPoint candidates[16];
    int found = 0;
    for (int i = 0; i < count; ++i) {
        // How far the incident shape's surface under this point is above the reference surface (both are
        // their radius out from the core). Negative: it has sunk in.
        const Real gap = dot(ref.normal, src[i].p - ref.center) - ref_radius - inc_radius;
        if (gap <= 0)
            candidates[found++] = {src[i].p - ref.normal * (inc_radius + gap * static_cast<Real>(0.5)), -gap, pair_code | (src[i].tag & 0x1FFFFFu)};
    }
    found = reduce_to_four(candidates, found, ref.u, ref.v);
    for (int i = 0; i < found; ++i) out[i] = candidates[i];
    return found;
}

}  // namespace

bool collide_convex(const Convex& a, const Convex& b, Manifold& out) {
    const ClosestResult hit = closest(a, b);
    if (hit.distance > 0) return false;

    // GJK and EPA give the direction and depth of the contact and one point of it. That is the whole
    // answer when the shapes meet at a point (a corner, two crossed edges, anything round).
    Manifold m;
    m.normal = hit.normal;
    m.depth = -hit.distance;
    m.count = 1;
    m.points[0] = {(hit.point_a + hit.point_b) * static_cast<Real>(0.5), m.depth, 0x80000000u};

    // But when a flat face is what is touching, the contact is a patch, and a body resting on one point
    // would rock. Look at the feature each shape presents to the other along the normal.
    const Feature fa = facing_feature(a, hit.normal), fb = facing_feature(b, -hit.normal);
    const bool face_a = fa.count == 4 && fa.alignment >= kFaceAligned, face_b = fb.count == 4 && fb.alignment >= kFaceAligned;
    if (face_a || face_b) {
        // If both present a face, a's is the reference unless b's is the squarer one by more than rounding
        // (so that two parallel faces do not swap roles from frame to frame).
        const bool ref_is_a = face_a && (!face_b || fa.alignment >= fb.alignment - static_cast<Real>(1e-6));
        ContactPoint patch[4];
        const int found = ref_is_a ? face_patch(fa, a.radius, fb, b.radius, true, patch) : face_patch(fb, b.radius, fa, a.radius, false, patch);
        if (found > 0) {
            // The face's own normal, exactly, in place of EPA's (they agree to within the 0.8 degrees
            // above). The depth has to go with the normal: how far the two overlap along THIS direction,
            // which the support functions give directly. (It can exceed every point's depth: the deepest
            // corner of the other shape may hang outside the face's outline.)
            m.normal = ref_is_a ? fa.normal : -fb.normal;
            m.depth = dot(m.normal, a.support(m.normal) - b.support(-m.normal));
            m.count = found;
            for (int i = 0; i < found; ++i) m.points[i] = patch[i];
        }
    } else if (fa.count == 2 && fb.count == 2) {
        // Two capsules side by side: if their axes are (nearly) parallel the contact is the stretch where
        // they run alongside each other, and its two ends are the contact points. Each end is measured for
        // itself, from a's axis straight across to b's, since "nearly" parallel axes are not the same
        // distance apart at both ends.
        const Vec3 mid_a = (fa.p[0] + fa.p[1]) * static_cast<Real>(0.5), along = (fa.p[1] - fa.p[0]).normalized();
        const Vec3 b0 = fb.p[0], b_axis = fb.p[1] - fb.p[0];
        const Real half_a = distance(fa.p[0], fa.p[1]) * static_cast<Real>(0.5);
        if (cross(along, b_axis.normalized()).length() < static_cast<Real>(0.02) && half_a > 0) {
            const Real t0 = dot(fb.p[0] - mid_a, along), t1 = dot(fb.p[1] - mid_a, along);
            const Real lo = std::max(std::min(t0, t1), -half_a), hi = std::min(std::max(t0, t1), half_a);
            ContactPoint ends[2];
            int found = 0;
            for (const Real t : {lo, hi}) {
                const Vec3 on_a = mid_a + along * t;
                const Vec3 on_b = b0 + b_axis * (dot(on_a - b0, b_axis) / b_axis.length_sq());
                const Real apart = distance(on_a, on_b), depth = a.radius + b.radius - apart;
                if (depth < 0 || apart == 0) continue;
                const Vec3 across = (on_b - on_a) / apart;
                ends[found] = {on_a + across * (a.radius - depth * static_cast<Real>(0.5)), depth, 0x40000000u + static_cast<std::uint32_t>(found)};
                ++found;
            }
            if (found == 2 && hi - lo > static_cast<Real>(1e-3) * half_a) {
                m.count = 2;
                m.points[0] = ends[0];
                m.points[1] = ends[1];
            }
        }
    }
    out = m;
    return true;
}

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
