#include <phys3d/gjk.hpp>

#include "vec3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace phys3d {

Convex Convex::sphere(Real r, Vec3 p) {
    Convex c;
    c.core = Core::Point;
    c.pos = p;
    c.radius = r;
    return c;
}
Convex Convex::capsule(Real half_length, Real r, Vec3 p, Quat orientation) {
    Convex c;
    c.core = Core::Segment;
    c.pos = p;
    c.q = orientation;
    c.half = {0, half_length, 0};
    c.radius = r;
    return c;
}
Convex Convex::box(Vec3 h, Vec3 p, Quat orientation) {
    Convex c;
    c.core = Core::Box;
    c.pos = p;
    c.q = orientation;
    c.half = h;
    return c;
}
Convex Convex::hull(const Vec3* pts, int n, Vec3 p, Quat orientation, Real r) {
    Convex c;
    c.core = Core::Points;
    c.pos = p;
    c.q = orientation;
    c.points = pts;
    c.count = n;
    c.radius = r;
    return c;
}
Convex Convex::hull(const Hull& built, Vec3 p, Quat orientation, Real r) {
    Convex c = hull(built.vertices.data(), static_cast<int>(built.vertices.size()), p, orientation, r);
    c.faces = &built;
    return c;
}
Convex Convex::of(const Body& body) {
    switch (body.shape.type) {
        case Shape::Type::Sphere: return sphere(body.shape.radius, body.pos);
        case Shape::Type::Capsule: return capsule(body.shape.half_length, body.shape.radius, body.pos, body.q);
        case Shape::Type::Hull: return hull(*body.shape.hull, body.pos, body.q);
        default: return box(body.shape.half_extents, body.pos, body.q);
    }
}

Vec3 Convex::core_corner(Vec3 d) const {
    if (core == Core::Point) return {};
    const Vec3 l = inv_rotate(q, d);  // the question, asked in the shape's own frame
    Vec3 best;
    switch (core) {
        case Core::Segment: best = {0, l.y < 0 ? -half.y : half.y, 0}; break;
        case Core::Box: best = {l.x < 0 ? -half.x : half.x, l.y < 0 ? -half.y : half.y, l.z < 0 ? -half.z : half.z}; break;
        default: {
            Real reach = -std::numeric_limits<Real>::max();
            for (int i = 0; i < count; ++i)
                if (dot(points[i], l) > reach) reach = dot(points[i], l), best = points[i];
        }
    }
    return best;
}

namespace {

// ---- In double precision ------------------------------------------------------------------------
//
// GJK and EPA work in double precision (step 18). Their answers are differences of support points, and
// next to a big shape those are big numbers: against a 300 m floor, single precision placed them to about
// 30 micrometres, so its tolerances had to be relative to the floor (a ten-thousandth of it, 4 cm) and a
// small shape near the floor read as touching it from 4 cm away. In double precision the same points are
// good to 10^-16 of their size, and the tolerances can be relative to the answer instead.
//
// That includes placing the support points. The corner is chosen in single precision, but turned and moved
// into place in double. In single precision a corner is placed to about a part in ten million of the shape's
// size (and, once the shape's position is added, of its distance from the world's origin: a millimetre at
// 10 km). That is far below anything the solver notices, but it decides which side of zero a contact that
// is exactly touching lands on. A crate resting on the floor tilted by 10^-7 radians has its lowest corner
// 9 * 10^-8 m into it; rounded to single precision and then measured in double, it came out 1.5 * 10^-8 m
// above, the contact began a step later, and a tower that had stood exactly as on the box code (step 10)
// ended 2 mm to one side of it. In double, GJK gets the 9 * 10^-8. (Choosing the corner in double too
// made no difference anyone could measure: a wrong choice is between two corners within 10^-7 of a tie.)

// The farthest point of the core in direction d, measured from the shape's position.
Vec3d core_offset(const Convex& c, Vec3d d) { return rotate_d(c.q, Vec3d(c.core_corner(to_real(d)))); }

// One sample of the Minkowski difference of the two cores, remembering which point of each it came from
// so the answer can be turned back into points on the shapes.
struct Vertex {
    Vec3d w, a, b;  // w = a - b
    Vec3 d;         // the direction it was sampled in (unit), kept for GjkCache
};

// Everything is done relative to `origin` (a's position), so the numbers stay small wherever in the
// world the pair happens to be. `shift` moves b, in double precision.
struct Pair {
    const Convex& a;
    const Convex& b;
    Vec3d origin;
    Vec3d shift{};
    Vec3d toward_b() const { return Vec3d(b.pos) + shift - Vec3d(a.pos); }
    Vertex sample(Vec3d d) const {
        Vertex v;
        v.a = (Vec3d(a.pos) - origin) + core_offset(a, d);
        v.b = (Vec3d(b.pos) - origin) + shift + core_offset(b, -d);
        v.w = v.a - v.b;
        v.d = to_real(d.normalized());
        return v;
    }
};

struct Simplex {
    Vertex v[4];
    double weight[4] = {1, 0, 0, 0};  // the closest point to the origin is sum(weight * v.w)
    int n = 0;

    Vec3d closest() const {
        Vec3d p;
        for (int i = 0; i < n; ++i) p += v[i].w * weight[i];
        return p;
    }
    double scale_sq() const {  // the size of the numbers involved, for tolerances
        double s = 0;
        for (int i = 0; i < n; ++i) s = std::max(s, v[i].w.length_sq());
        return s;
    }
};

void closest_on_segment(Vec3d a, Vec3d b, double out[2]) {
    const Vec3d ab = b - a;
    const double len_sq = ab.length_sq();
    const double t = len_sq > 0 ? std::clamp(dot(-a, ab) / len_sq, 0.0, 1.0) : 0.0;
    out[0] = 1 - t;
    out[1] = t;
}

// The plane around a triangle splits into seven regions (three corners, three edges, the inside); each
// test below asks "is the origin in this one?" using only dot products. (Ericson, Real-Time Collision
// Detection, 5.1.5.) The tests multiply dot products together, so their size goes as the FOURTH power of
// the triangle's: for the samples of a 60 m floor against a ball 20 cm above it that is 10^7, which is why
// this routine was the first part of GJK to move to double precision (step 17).
void closest_on_triangle(Vec3d a, Vec3d b, Vec3d c, double out[3]) {
    const Vec3d ab = b - a, ac = c - a;
    const double d1 = dot(ab, -a), d2 = dot(ac, -a);
    out[0] = out[1] = out[2] = 0;
    if (d1 <= 0 && d2 <= 0) { out[0] = 1; return; }  // corner a
    const double d3 = dot(ab, -b), d4 = dot(ac, -b);
    if (d3 >= 0 && d4 <= d3) { out[1] = 1; return; }  // corner b
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {  // edge ab
        const double t = d1 / (d1 - d3);
        out[0] = 1 - t, out[1] = t;
        return;
    }
    const double d5 = dot(ab, -c), d6 = dot(ac, -c);
    if (d6 >= 0 && d5 <= d6) { out[2] = 1; return; }  // corner c
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {  // edge ac
        const double t = d2 / (d2 - d6);
        out[0] = 1 - t, out[2] = t;
        return;
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) {  // edge bc
        const double t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        out[1] = 1 - t, out[2] = t;
        return;
    }
    const double sum = va + vb + vc;
    out[1] = vb / sum;
    out[2] = vc / sum;
    out[0] = 1 - out[1] - out[2];
}

}  // namespace

namespace detail {

void closest_on_segment(Vec3 a, Vec3 b, Real out[2]) {
    double w[2];
    phys3d::closest_on_segment(Vec3d(a), Vec3d(b), w);
    for (int i = 0; i < 2; ++i) out[i] = static_cast<Real>(w[i]);
}

void closest_on_triangle(Vec3 a, Vec3 b, Vec3 c, Real out[3]) {
    double w[3];
    phys3d::closest_on_triangle(Vec3d(a), Vec3d(b), Vec3d(c), w);
    for (int i = 0; i < 3; ++i) out[i] = static_cast<Real>(w[i]);
}

}  // namespace detail

namespace {

// Finds the point of the simplex closest to the origin, sets the weights, and throws away the vertices
// that do not support it (so a tetrahedron whose nearest feature is an edge becomes that edge). Returns
// true if the simplex is a tetrahedron with the origin inside.
bool solve(Simplex& s) {
    double w[4] = {0, 0, 0, 0};
    if (s.n == 1) {
        w[0] = 1;
    } else if (s.n == 2) {
        closest_on_segment(s.v[0].w, s.v[1].w, w);
    } else if (s.n == 3) {
        closest_on_triangle(s.v[0].w, s.v[1].w, s.v[2].w, w);
    } else {
        // The origin is inside unless it is on the outer side of some face, the side away from the fourth
        // vertex. For each such face, the nearest point of that face is a candidate.
        static constexpr int kFace[4][4] = {{0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {1, 3, 2, 0}};  // three corners, then the opposite vertex
        double best = std::numeric_limits<double>::max();
        bool outside_any = false;
        for (const auto& f : kFace) {
            const Vec3d a = s.v[f[0]].w, b = s.v[f[1]].w, c = s.v[f[2]].w, opposite = s.v[f[3]].w;
            const Vec3d n = cross(b - a, c - a);
            const double origin_side = dot(n, -a), inner_side = dot(n, opposite - a);
            // A flat tetrahedron encloses nothing, but rounding can make each of these sign tests pass for
            // one, and GJK would then announce an overlap between shapes that are a hair apart. So a
            // tetrahedron whose fourth vertex is within a billionth of its size of the opposite face counts
            // as flat, and every face of it as outer. (inner_side / |n| is that vertex's height.)
            const double flat = 1e-9 * n.length() * std::sqrt(std::max({(b - a).length_sq(), (c - a).length_sq(), (opposite - a).length_sq()}));
            if (std::fabs(inner_side) > flat && origin_side * inner_side >= 0) continue;
            outside_any = true;
            double t[3];
            closest_on_triangle(a, b, c, t);
            const double dist_sq = (a * t[0] + b * t[1] + c * t[2]).length_sq();
            if (dist_sq < best) {
                best = dist_sq;
                w[0] = w[1] = w[2] = w[3] = 0;
                w[f[0]] = t[0], w[f[1]] = t[1], w[f[2]] = t[2];
            }
        }
        if (!outside_any) return true;
    }
    int kept = 0;
    for (int i = 0; i < s.n; ++i)
        if (w[i] > 0) {
            s.v[kept] = s.v[i];
            s.weight[kept] = w[i];
            ++kept;
        }
    s.n = kept;
    return false;
}

// The unit direction from the origin to the simplex's nearest point, for a simplex that does not contain
// the origin. Adding up weight * vertex gives that point, but when the shapes are a hair apart the point
// is a tiny difference of large numbers and its direction is mostly rounding. For a triangle the geometry
// gives the direction without that cancellation: the nearest point of its interior lies along its normal.
Vec3d direction_to(const Simplex& s) {
    if (s.n != 3) return s.closest().normalized();
    const Vec3d n = cross(s.v[1].w - s.v[0].w, s.v[2].w - s.v[0].w).normalized();
    return dot(n, s.v[0].w) < 0 ? -n : n;
}

constexpr int kMaxGjkIterations = 48;

// Runs GJK on the cores. Returns true if they overlap (the simplex then holds or encloses the origin);
// otherwise the simplex and its weights give the closest points.
bool gjk(const Pair& pair, Simplex& s, int& iterations, bool& converged, const GjkCache* cache = nullptr) {
    converged = true;
    iterations = 1;
    Vec3d v;
    if (cache && cache->count > 0) {
        // Start from last step's simplex: sample A - B again in the directions that found its corners.
        // The shapes have barely moved, so these are nearly the corners GJK would end at, and for a pair
        // that was touching they usually still enclose the origin, which is the whole of GJK's job.
        s.n = 0;
        for (int k = 0; k < cache->count; ++k) {
            const Vertex x = pair.sample(Vec3d(cache->directions[k]));
            bool repeat = false;
            for (int i = 0; i < s.n; ++i) repeat = repeat || (s.v[i].w - x.w).length_sq() == 0;
            if (!repeat) s.v[s.n++] = x;
        }
        if (solve(s)) return true;
        v = s.closest();
    } else {
        // Otherwise look first toward b from a.
        Vec3d start = pair.toward_b();
        if (start.length_sq() == 0) start = {1, 0, 0};
        s.v[0] = pair.sample(start);
        s.weight[0] = 1;
        s.n = 1;
        v = s.v[0].w;
    }

    for (iterations = 1; iterations <= kMaxGjkIterations; ++iterations) {
        const double v_sq = v.length_sq();
        const double scale_sq = s.scale_sq();
        // Closer than a ten-billionth of the shapes' size counts as touching, which counts as overlap: the
        // direction of a v that short is mostly rounding, and EPA finds the normal instead. (In single
        // precision this had to be a ten-thousandth.)
        if (v_sq <= 1e-20 * scale_sq) return true;

        const Vertex w = pair.sample(-v);
        // v is the nearest point found so far, at distance |v|. The plane through w perpendicular to v has
        // all of A - B on its far side, so the true distance is at least (v . w) / |v|. When the two agree
        // we are done. (Agreeing to one part in a million in distance squared pins the distance, but the
        // direction of v only to about a thousandth of a radian. The cores are polytopes, though, and GJK
        // nearly always ends on the exact answer first, by sampling a corner it already has: tightening
        // this to 10^-10 changed nothing measurable.)
        if (v_sq - dot(v, w.w) <= 1e-6 * v_sq) return false;
        for (int i = 0; i < s.n; ++i)
            if ((s.v[i].w - w.w).length_sq() == 0) return false;  // the same point again: nothing further that way

        const Simplex before = s;
        s.v[s.n] = w;
        ++s.n;
        if (solve(s)) return true;
        const Vec3d next = s.closest();
        if (!(next.length_sq() < v_sq)) {  // rounding: no nearer than before. Keep the better answer.
            s = before;
            return false;
        }
        v = next;
    }
    converged = false;
    return false;
}

// ---- EPA -----------------------------------------------------------------------------------------

constexpr int kMaxEpaVertices = 80;
constexpr int kMaxEpaFaces = 2 * kMaxEpaVertices;  // a closed triangle mesh with V vertices has 2V - 4 faces
constexpr int kMaxEpaIterations = kMaxEpaVertices - 4;

struct Face {
    int i[3];
    Vec3d n;      // unit, pointing away from the origin
    double dist;  // of the face's plane from the origin
};

struct Polytope {
    Vertex v[kMaxEpaVertices];
    int vertex_count = 0;
    Face f[kMaxEpaFaces];
    int face_count = 0;

    // Adds the triangle (i0, i1, i2), wound counter-clockwise seen from outside. False if it has no area.
    bool add_face(int i0, int i1, int i2) {
        const Vec3d a = v[i0].w, e1 = v[i1].w - a, e2 = v[i2].w - a;
        const Vec3d n = cross(e1, e2);
        const double len = n.length();
        if (!(len > 0) || face_count == kMaxEpaFaces) return false;
        f[face_count++] = {{i0, i1, i2}, n / len, dot(n, a) / len};
        return true;
    }
};

// Distance of p from the line through a and b.
double distance_from_line(Vec3d p, Vec3d a, Vec3d b) {
    const Vec3d ab = b - a;
    return cross(p - a, ab).length() / std::max(ab.length(), 1e-300);
}

// Some direction perpendicular to the unit vector `along`.
Vec3d perpendicular(Vec3d along) { return cross(along, std::fabs(along.x) < 0.57735 ? Vec3d{1, 0, 0} : Vec3d{0, 1, 0}).normalized(); }

struct Penetration {
    Vec3d normal;  // from a toward b
    double depth = 0;
    int iterations = 0;
    bool converged = true;
};

// The cores overlap and `s` holds or encloses the origin. Finds the shortest way out of A - B.
Penetration epa(const Pair& pair, Simplex s) {
    Penetration out;
    // How big A - B is, for tolerances. The simplex's own points say, unless they all sit on the origin
    // (the cores touch at one point), in which case sample along the axes.
    double size_sq = s.scale_sq();
    if (s.n < 4)
        for (const Vec3d& d : {Vec3d{1, 0, 0}, Vec3d{0, 1, 0}, Vec3d{0, 0, 1}, Vec3d{-1, 0, 0}, Vec3d{0, -1, 0}, Vec3d{0, 0, -1}})
            size_sq = std::max(size_sq, pair.sample(d).w.length_sq());
    const double size = std::sqrt(size_sq);
    const double apart = 1e-4 * size;  // two samples closer than this are "the same"

    // EPA needs a real tetrahedron to start from, and GJK may have stopped with less: the origin lay
    // exactly on a vertex, edge or face. Build it up by sampling in directions the simplex does not span.
    if (s.n == 1) {
        const Vec3d axes[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
        for (const Vec3d& d : axes) {
            const Vertex c = pair.sample(d);
            if (distance(c.w, s.v[0].w) > apart) {
                s.v[s.n++] = c;
                break;
            }
        }
    }
    if (s.n == 2) {
        // Six directions around the segment, 60 degrees apart.
        const Vec3d along = (s.v[1].w - s.v[0].w).normalized();
        const Vec3d side = perpendicular(along), third = cross(along, side);
        for (int k = 0; k < 6; ++k) {
            const double angle = k * (3.14159265358979323846 / 3);
            const Vertex c = pair.sample(side * std::cos(angle) + third * std::sin(angle));
            if (distance_from_line(c.w, s.v[0].w, s.v[1].w) > apart) {
                s.v[s.n++] = c;
                break;
            }
        }
    }
    if (s.n == 3) {
        const Vec3d n = cross(s.v[1].w - s.v[0].w, s.v[2].w - s.v[0].w).normalized();
        for (const Vec3d& d : {n, -n}) {
            const Vertex c = pair.sample(d);
            if (std::fabs(dot(n, c.w - s.v[0].w)) > apart) {
                s.v[s.n++] = c;
                break;
            }
        }
    }
    if (s.n < 4) {
        // A - B is flat (or a line, or a point): two segments crossing, a point on a face. The cores touch
        // without any depth, and the way apart is perpendicular to whatever A - B spans.
        const Vec3d toward_b = pair.toward_b();
        Vec3d n;
        if (s.n == 3) {
            n = cross(s.v[1].w - s.v[0].w, s.v[2].w - s.v[0].w).normalized();
        } else if (s.n == 2) {
            const Vec3d along = (s.v[1].w - s.v[0].w).normalized();
            n = (toward_b - along * dot(toward_b, along)).normalized();
            if (n.length_sq() == 0) n = perpendicular(along);
        } else {
            n = toward_b.normalized();
            if (n.length_sq() == 0) n = {0, 1, 0};
        }
        out.normal = n;
        return out;
    }

    Polytope poly;
    for (int i = 0; i < 4; ++i) poly.v[i] = s.v[i];
    poly.vertex_count = 4;
    // Make vertex 3 lie behind face (0, 1, 2), so the four faces below are all wound outward.
    if (dot(cross(poly.v[1].w - poly.v[0].w, poly.v[2].w - poly.v[0].w), poly.v[3].w - poly.v[0].w) > 0) std::swap(poly.v[0], poly.v[1]);
    bool ok = poly.add_face(0, 1, 2) && poly.add_face(0, 3, 1) && poly.add_face(0, 2, 3) && poly.add_face(1, 3, 2);

    Face nearest{};
    auto find_nearest = [&] {
        nearest = poly.f[0];
        for (int k = 1; k < poly.face_count; ++k)
            if (poly.f[k].dist < nearest.dist) nearest = poly.f[k];
    };
    if (ok) find_nearest();
    out.converged = false;
    for (out.iterations = 1; ok && out.iterations <= kMaxEpaIterations; ++out.iterations) {
        // Does A - B reach beyond the nearest face? If not, that face is on its surface.
        const Vertex w = pair.sample(nearest.n);
        if (dot(w.w, nearest.n) - nearest.dist <= 1e-9 * size) {
            out.converged = true;
            break;
        }
        if (poly.vertex_count == kMaxEpaVertices) break;
        const int added = poly.vertex_count;
        poly.v[poly.vertex_count++] = w;

        // Remove every face the new point can see. What is left has a hole whose rim (the horizon) is the
        // edges that belonged to exactly one removed face: an edge shared by two removed faces shows up
        // once in each direction and cancels.
        int edge[kMaxEpaFaces * 3][2];
        int edges = 0;
        for (int k = 0; k < poly.face_count;) {
            const Face& face = poly.f[k];
            if (dot(face.n, w.w) - face.dist <= 0) {
                ++k;
                continue;
            }
            for (int e = 0; e < 3; ++e) {
                const int from = face.i[e], to = face.i[(e + 1) % 3];
                bool cancelled = false;
                for (int j = 0; j < edges && !cancelled; ++j)
                    if (edge[j][0] == to && edge[j][1] == from) {
                        edge[j][0] = edge[edges - 1][0];
                        edge[j][1] = edge[edges - 1][1];
                        --edges;
                        cancelled = true;
                    }
                if (!cancelled) edge[edges][0] = from, edge[edges][1] = to, ++edges;
            }
            poly.f[k] = poly.f[--poly.face_count];
        }
        // Close the hole with a fan of new faces from the rim to the new point.
        for (int j = 0; j < edges && ok; ++j) ok = poly.add_face(edge[j][0], edge[j][1], added);
        if (ok) find_nearest();  // otherwise keep the last good answer
    }

    out.normal = nearest.n;
    out.depth = std::max(nearest.dist, 0.0);
    return out;
}

// From the cores to the rounded shapes: each surface is its radius further out along the normal.
void place_points(const Convex& a, const Convex& b, const Simplex& s, Vec3d shift, ClosestResult& out) {
    Vec3d core_a, core_b;
    for (int i = 0; i < s.n; ++i) {
        core_a += s.v[i].a * s.weight[i];
        core_b += s.v[i].b * s.weight[i];
    }
    const Vec3d origin(a.pos);
    out.point_a = to_real(origin + core_a) + out.normal * a.radius;
    out.point_b = to_real(origin + core_b - shift) - out.normal * b.radius;
    out.has_points = true;
}

}  // namespace

void find_deepest_points(const Convex& a, const Convex& b, ClosestResult& out) {
    if (out.has_points) return;
    // Where are the deepest points? EPA's nearest face could tell us, but only when the face it stopped on
    // happens to contain the foot of the perpendicular from the origin, and a flat side of A - B is covered
    // by several triangles (or, right at the start, by part of one). So ask the question GJK answers
    // reliably instead: pull b out along the normal until the cores are just clear, and take the closest
    // points of that pair. "Just clear" is a thousandth of the shapes' size.
    const double core_depth = -(static_cast<double>(out.distance) + a.radius + b.radius);
    const double clearance = 1e-3 * std::max(static_cast<double>(out.size), core_depth);
    // (b is moved in double precision: moving the Convex itself would round its new position to single
    // precision, by up to 30 micrometres at 300 m from the origin.)
    const Vec3d shift = Vec3d(out.normal) * (core_depth + clearance);
    Simplex s;
    int extra = 0;
    bool extra_converged = true;
    // The pulled-apart pair should be exactly `clearance` apart. If GJK says otherwise (or says they still
    // overlap), one of the two runs stopped early on rounding. The depth and normal still hold, but the
    // points may be out, so say so.
    if (gjk(Pair{a, b, Vec3d(a.pos), shift}, s, extra, extra_converged) || std::fabs(dot(direction_to(s), s.v[0].w) - clearance) > 0.5 * clearance)
        out.converged = false;
    out.converged = out.converged && extra_converged;
    out.gjk_iterations += extra;
    place_points(a, b, s, shift, out);
}
ClosestResult closest(const Convex& a, const Convex& b, bool with_points, GjkCache* cache) {
    const Pair pair{a, b, Vec3d(a.pos)};
    ClosestResult out;
    Simplex s;
    const bool cores_overlap = gjk(pair, s, out.gjk_iterations, out.converged, cache);
    if (cache) {
        cache->count = s.n;
        cache->iterations = out.gjk_iterations;
        for (int i = 0; i < s.n; ++i) cache->directions[i] = s.v[i].d;
    }
    if (cores_overlap) {
        // How big the shapes are, for the tolerance used if the deepest points are asked for: the farthest
        // of the points of either that GJK sampled, measured from a's position. (Not the samples of A - B
        // themselves, which can all sit at the origin when the cores only touch.)
        double size_sq = 0;
        for (int i = 0; i < s.n; ++i) size_sq = std::max({size_sq, s.v[i].a.length_sq(), s.v[i].b.length_sq()});
        out.size = static_cast<Real>(std::sqrt(size_sq));
        const Penetration p = epa(pair, s);
        out.epa_iterations = p.iterations;
        out.converged = out.converged && p.converged;
        out.normal = to_real(p.normal);
        out.distance = static_cast<Real>(-p.depth - a.radius - b.radius);
        out.has_points = false;
        if (with_points) find_deepest_points(a, b, out);
        return out;
    }
    // Apart. The nearest point of A - B to the origin is (a point of a) - (a point of b), so it points from
    // b back to a, and its length is the distance.
    const Vec3d toward_a = direction_to(s);
    out.normal = to_real(-toward_a);
    out.distance = static_cast<Real>(dot(toward_a, s.v[0].w) - a.radius - b.radius);
    place_points(a, b, s, {}, out);
    return out;
}

}  // namespace phys3d
