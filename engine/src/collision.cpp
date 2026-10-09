#include <phys/collision.hpp>

#include <algorithm>
#include <limits>

namespace phys {

namespace {

constexpr Real kRealMax = std::numeric_limits<Real>::max();

// When two faces have almost the same separation, keep using A's face as the reference. Without this
// the choice flips between A and B from frame to frame and contact points jitter.
constexpr Real kReferenceFaceBias = static_cast<Real>(0.0005);

// ---------------------------------------------------------------------------------------------
// circle vs circle
// ---------------------------------------------------------------------------------------------
bool collide_circles(const Body& a, const Body& b, Manifold& out) {
    const Real ra = a.shape.circle.radius, rb = b.shape.circle.radius;
    const Vec2 d = b.pos - a.pos;
    const Real dist_sq = d.length_sq();
    const Real r = ra + rb;
    if (dist_sq > r * r) return false;

    const Real dist = std::sqrt(dist_sq);
    // Concentric circles have no preferred direction; any unit vector is a valid way out.
    const Vec2 n = dist > static_cast<Real>(1e-9) ? d * (1 / dist) : Vec2{1, 0};
    out.normal = n;
    out.depth = r - dist;
    out.count = 1;
    out.points[0] = {a.pos + n * (ra - out.depth * static_cast<Real>(0.5)), out.depth};
    return true;
}

// ---------------------------------------------------------------------------------------------
// polygon vs circle (normal points from the polygon to the circle)
// ---------------------------------------------------------------------------------------------
bool collide_polygon_circle(const Body& poly_body, const Body& circle_body, Manifold& out) {
    const Polygon& p = poly_body.shape.polygon;
    const Real r = circle_body.shape.circle.radius;

    // Work in the polygon's own frame, where it has no rotation.
    const Vec2 c = apply_inv(poly_body.transform(), circle_body.pos);

    // Face of greatest separation: how far the circle centre is in front of each face plane.
    int face = 0;
    Real sep = -kRealMax;
    for (int i = 0; i < p.count; ++i) {
        const Real s = dot(p.normals[i], c - p.vertices[i]);
        if (s > r) return false;  // the circle is entirely in front of this face
        if (s > sep) {
            sep = s;
            face = i;
        }
    }

    const Vec2 v1 = p.vertices[face], v2 = p.vertices[(face + 1) % p.count];
    Vec2 n_local;
    Real s;  // signed distance from the polygon surface to the circle centre, along n_local
    std::uint32_t id;  // the polygon feature touched: a face index, or kVertexFlag | vertex index
    constexpr std::uint32_t kVertexFlag = 0x100;
    if (sep < std::numeric_limits<Real>::epsilon()) {
        // Centre is inside the polygon: push out through the nearest face.
        n_local = p.normals[face];
        s = sep;
        id = static_cast<std::uint32_t>(face);
    } else if (dot(c - v1, v2 - v1) <= 0) {
        // Nearest feature is vertex v1.
        const Vec2 d = c - v1;
        s = d.length();
        if (s > r) return false;
        n_local = d.normalized();
        id = kVertexFlag | static_cast<std::uint32_t>(face);
    } else if (dot(c - v2, v1 - v2) <= 0) {
        // Nearest feature is vertex v2.
        const Vec2 d = c - v2;
        s = d.length();
        if (s > r) return false;
        n_local = d.normalized();
        id = kVertexFlag | static_cast<std::uint32_t>((face + 1) % p.count);
    } else {
        // Nearest feature is the face itself.
        n_local = p.normals[face];
        s = sep;
        id = static_cast<std::uint32_t>(face);
    }

    const Vec2 n = rotate(poly_body.q, n_local);
    out.normal = n;
    out.depth = r - s;
    out.count = 1;
    // Circle's deepest point is centre - n*r, the polygon surface is centre - n*s; take the midpoint.
    out.points[0] = {circle_body.pos - n * ((r + s) * static_cast<Real>(0.5)), out.depth, id};
    return true;
}

// ---------------------------------------------------------------------------------------------
// polygon vs polygon: separating axis test, then clip to build the manifold
// ---------------------------------------------------------------------------------------------
struct WorldPolygon {
    int count = 0;
    Vec2 v[Polygon::kMaxVertices];
    Vec2 n[Polygon::kMaxVertices];
};

WorldPolygon to_world(const Body& b) {
    const Polygon& p = b.shape.polygon;
    WorldPolygon w;
    w.count = p.count;
    for (int i = 0; i < p.count; ++i) {
        w.v[i] = apply(b.transform(), p.vertices[i]);
        w.n[i] = rotate(b.q, p.normals[i]);
    }
    return w;
}

// For each face of `a`, measure how far the nearest vertex of `b` is in front of that face's plane
// (negative = behind it). The face with the LARGEST such value is the axis of least penetration.
// If it is positive, that face is a separating axis and the polygons do not touch.
Real find_max_separation(const WorldPolygon& a, const WorldPolygon& b, int& best_face) {
    Real best = -kRealMax;
    best_face = 0;
    for (int i = 0; i < a.count; ++i) {
        Real min_sep = kRealMax;
        for (int j = 0; j < b.count; ++j) min_sep = std::min(min_sep, dot(a.n[i], b.v[j] - a.v[i]));
        if (min_sep > best) {
            best = min_sep;
            best_face = i;
        }
    }
    return best;
}

// A point on the incident edge plus a tag naming where it came from: 0..7 is a vertex of the
// incident polygon, kClipTag + ... is a point created by clipping against a side plane. The tag is
// what lets us recognise "the same" contact point in the next frame.
struct ClipVertex {
    Vec2 p;
    std::uint8_t tag = 0;
};
constexpr std::uint8_t kClipTag = 16;
constexpr std::uint8_t kFallbackTag = 40;

// Sutherland-Hodgman against one half-plane: keep the part of segment `in` where dot(n, p) <= c.
// `plane` (0 or 1) says which side plane this is, for tagging. Returns how many points survive.
int clip_segment(ClipVertex out[2], const ClipVertex in[2], Vec2 n, Real c, int plane) {
    int count = 0;
    const Real d0 = dot(n, in[0].p) - c, d1 = dot(n, in[1].p) - c;
    if (d0 <= 0) out[count++] = in[0];
    if (d1 <= 0) out[count++] = in[1];
    if (d0 * d1 < 0) {  // endpoints on opposite sides: add the crossing point
        const Real alpha = d0 / (d0 - d1);
        ClipVertex v;
        v.p = in[0].p + (in[1].p - in[0].p) * alpha;
        // The crossing replaces whichever endpoint was cut off.
        v.tag = static_cast<std::uint8_t>(kClipTag + 2 * plane + (d0 > 0 ? 0 : 1));
        out[count++] = v;
    }
    return count;
}

bool collide_polygons(const Body& a_body, const Body& b_body, Manifold& out) {
    const WorldPolygon A = to_world(a_body), B = to_world(b_body);

    int face_a, face_b;
    const Real sep_a = find_max_separation(A, B, face_a);
    if (sep_a > 0) return false;
    const Real sep_b = find_max_separation(B, A, face_b);
    if (sep_b > 0) return false;

    // The reference polygon supplies the contact normal (its least-penetrated face); the other,
    // the incident polygon, supplies the points that poke through.
    const WorldPolygon *ref = &A, *inc = &B;
    int ref_face = face_a;
    Real sep = sep_a;
    bool flipped = false;
    if (sep_b > sep_a + kReferenceFaceBias) {
        ref = &B;
        inc = &A;
        ref_face = face_b;
        sep = sep_b;
        flipped = true;
    }
    const Vec2 ref_n = ref->n[ref_face];

    // Incident edge: the edge of the other polygon facing most directly against the reference face.
    int inc_face = 0;
    Real min_d = kRealMax;
    for (int j = 0; j < inc->count; ++j) {
        const Real d = dot(ref_n, inc->n[j]);
        if (d < min_d) {
            min_d = d;
            inc_face = j;
        }
    }
    const ClipVertex incident[2] = {
        {inc->v[inc_face], static_cast<std::uint8_t>(inc_face)},
        {inc->v[(inc_face + 1) % inc->count], static_cast<std::uint8_t>((inc_face + 1) % inc->count)}};

    // Clip the incident edge to the "side planes" of the reference edge, so only the part that lies
    // over the reference face can produce contact points.
    const Vec2 v1 = ref->v[ref_face], v2 = ref->v[(ref_face + 1) % ref->count];
    const Vec2 t = (v2 - v1).normalized();
    ClipVertex clipped_a[2], clipped_b[2];
    const bool clipped = clip_segment(clipped_a, incident, -t, -dot(t, v1), 0) == 2 &&
                         clip_segment(clipped_b, clipped_a, t, dot(t, v2), 1) == 2;

    // A point's id combines the two edges in contact with the tag of the point on the incident edge,
    // so it changes only when a different pair of features touches.
    auto make_id = [&](std::uint8_t tag) {
        return (flipped ? 1u << 24 : 0u) | (static_cast<std::uint32_t>(ref_face) << 16) |
               (static_cast<std::uint32_t>(inc_face) << 8) | tag;
    };

    out.normal = flipped ? -ref_n : ref_n;  // always A -> B
    out.depth = -sep;
    out.count = 0;
    if (clipped) {
        for (int k = 0; k < 2; ++k) {
            const Real s = dot(ref_n, clipped_b[k].p - v1);  // <= 0 means behind the reference face
            if (s <= 0)
                out.points[out.count++] = {clipped_b[k].p - ref_n * (s * static_cast<Real>(0.5)), -s,
                                           make_id(clipped_b[k].tag)};
        }
    }
    if (out.count == 0) {
        // The polygons overlap (SAT proved it) but the clipped edge never dips behind the reference
        // face; this happens when they just clip each other's corners. Report the deepest vertex of
        // the incident polygon so the pair is still pushed apart along the correct normal.
        const Real s0 = dot(ref_n, incident[0].p - v1), s1 = dot(ref_n, incident[1].p - v1);
        const Vec2 deepest = s1 < s0 ? incident[1].p : incident[0].p;
        const Real s = std::min(std::min(s0, s1), Real(0));
        out.points[out.count++] = {deepest - ref_n * (s * static_cast<Real>(0.5)), out.depth,
                                   make_id(kFallbackTag)};
    }
    return true;
}

}  // namespace

bool collide(const Body& a, const Body& b, Manifold& out) {
    using Type = Shape::Type;
    Manifold m;
    bool hit;
    if (a.shape.type == Type::Circle && b.shape.type == Type::Circle) {
        hit = collide_circles(a, b, m);
    } else if (a.shape.type == Type::Polygon && b.shape.type == Type::Circle) {
        hit = collide_polygon_circle(a, b, m);
    } else if (a.shape.type == Type::Circle && b.shape.type == Type::Polygon) {
        hit = collide_polygon_circle(b, a, m);
        m.normal = -m.normal;  // computed polygon -> circle; we promise A -> B
    } else {
        hit = collide_polygons(a, b, m);
    }
    if (hit) out = m;
    return hit;
}

}  // namespace phys
