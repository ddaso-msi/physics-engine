#include <phys3d/hull.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace phys3d {

namespace {

// The builder does its own arithmetic in double precision. Its tolerance is a hundred-thousandth of the
// cloud's size, and single precision cannot place the plane of a thin triangle that well: a triangle a
// hundredth as wide as it is long has a normal good to about one part in 100000, which is the whole
// tolerance gone. The engine's Real is used for what goes in and what comes out.
struct Vec {
    double x = 0, y = 0, z = 0;
    Vec() = default;
    Vec(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    explicit Vec(Vec3 v) : x(v.x), y(v.y), z(v.z) {}
    double length() const { return std::sqrt(x * x + y * y + z * z); }
};
Vec operator-(Vec a, Vec b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec operator*(Vec a, double s) { return {a.x * s, a.y * s, a.z * s}; }
double dot(Vec a, Vec b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec cross(Vec a, Vec b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
double distance(Vec a, Vec b) { return (a - b).length(); }

struct Triangle {
    int v[3];     // indices into the input points, counter-clockwise seen from outside
    Vec normal;   // unit
    double offset;
    bool alive = true;
};

// Distance of p from the line through a and b.
double distance_from_line(Vec p, Vec a, Vec b) {
    const Vec ab = b - a;
    return cross(p - a, ab).length() / ab.length();
}

struct Builder {
    std::vector<Vec> points;  // the input, in double precision
    int count = 0;
    double eps = 0;    // how far off a plane still counts as on it, for deciding what is a face
    double exact = 0;  // rounding in the builder's own arithmetic: anything under this is zero
    std::vector<Triangle> triangles;

    const Vec& at(int i) const { return points[static_cast<size_t>(i)]; }

    // Adds the triangle (a, b, c). False if it has no area.
    bool add(int a, int b, int c) {
        const Vec n = cross(at(b) - at(a), at(c) - at(a));
        const double len = n.length();
        if (!(len > 0)) return false;
        triangles.push_back({{a, b, c}, n * (1 / len), dot(n, at(a)) / len, true});
        return true;
    }
    double height(const Triangle& t, int p) const { return dot(t.normal, at(p)) - t.offset; }

    // Four points spread as widely as possible: the two ends of the x range, the point farthest from the
    // line through them, and the point farthest from the plane through those three.
    bool start() {
        int lo = 0, hi = 0;
        for (int i = 1; i < count; ++i) {
            if (at(i).x < at(lo).x) lo = i;
            if (at(i).x > at(hi).x) hi = i;
        }
        if (!(distance(at(lo), at(hi)) > eps)) return false;

        int third = -1;
        double best = eps;
        for (int i = 0; i < count; ++i) {
            const double d = distance_from_line(at(i), at(lo), at(hi));
            if (d > best) best = d, third = i;
        }
        if (third < 0) return false;

        Vec n = cross(at(hi) - at(lo), at(third) - at(lo));
        n = n * (1 / n.length());
        int fourth = -1;
        best = eps;
        for (int i = 0; i < count; ++i) {
            const double d = std::fabs(dot(n, at(i) - at(lo)));
            if (d > best) best = d, fourth = i;
        }
        if (fourth < 0) return false;

        // Wind the first triangle so that the fourth point is behind it; the other three follow.
        int a = lo, b = hi;
        if (dot(n, at(fourth) - at(lo)) > 0) std::swap(a, b);
        return add(a, b, third) && add(a, fourth, b) && add(a, third, fourth) && add(b, fourth, third);
    }

    // The live triangle on the other side of the edge from -> to (the one that has it as to -> from), or -1.
    int across(int from, int to) const {
        for (size_t i = 0; i < triangles.size(); ++i)
            if (triangles[i].alive)
                for (int k = 0; k < 3; ++k)
                    if (triangles[i].v[k] == to && triangles[i].v[(k + 1) % 3] == from) return static_cast<int>(i);
        return -1;
    }

    // Brings point p into the hull, given one triangle it is in front of. False if the surface could not be
    // closed around it.
    bool absorb(int p, int seen_from) {
        // The patch of triangles p sees: spread out from the first across shared edges, taking every
        // neighbour p is in front of. ("In front of" means by more than the arithmetic's own rounding. With
        // the points held in double precision that is as good as exact, and exactness is what keeps the
        // patch a disc with a single rim. An earlier version judged this in single precision with a
        // tolerance, and whichever way the tolerance leaned, something broke: lean toward "seen" and the
        // surface shrank away from corners it already had; lean toward "not seen" and new triangles were
        // laid on top of old ones.)
        std::vector<int> removed{seen_from}, frontier{seen_from};
        triangles[static_cast<size_t>(seen_from)].alive = false;
        while (!frontier.empty()) {
            const Triangle t = triangles[static_cast<size_t>(frontier.back())];
            frontier.pop_back();
            for (int e = 0; e < 3; ++e) {
                const int other = across(t.v[e], t.v[(e + 1) % 3]);
                if (other < 0 || !(height(triangles[static_cast<size_t>(other)], p) > exact)) continue;
                triangles[static_cast<size_t>(other)].alive = false;
                removed.push_back(other);
                frontier.push_back(other);
            }
        }
        // Its rim: the edges of removed triangles that still have a live triangle across them. The hole is
        // closed by a triangle from each rim edge to p.
        //
        // (Could p be exactly in line with a rim edge, making a triangle with no area? Not when points are
        // brought in farthest first. If p lay on the line through corners A and B beyond B, it would be at
        // least twice as far outside every triangle as B ever was, and would have been brought in before B.)
        std::vector<std::pair<int, int>> rim;
        for (const int t : removed)
            for (int e = 0; e < 3; ++e) {
                const int from = triangles[static_cast<size_t>(t)].v[e], to = triangles[static_cast<size_t>(t)].v[(e + 1) % 3];
                if (across(from, to) >= 0) rim.emplace_back(from, to);
            }
        for (const auto& edge : rim)
            if (!add(edge.first, edge.second, p)) return false;
        return true;
    }
};

bool build(const Vec3* points, int count, Hull& out) {
    if (count < 4) return false;

    Builder b;
    b.count = count;
    Vec3 lo = points[0], hi = points[0];
    for (int i = 0; i < count; ++i) {
        b.points.emplace_back(points[i]);
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], points[i][k]);
            hi[k] = std::max(hi[k], points[i][k]);
        }
    }
    b.eps = 1e-5 * static_cast<double>(distance(lo, hi));
    b.exact = 1e-10 * static_cast<double>(distance(lo, hi));
    if (!b.start()) return false;

    // Farthest point first, until nothing is outside. (Farthest first means the hull grows in big steps and
    // most points are swallowed without ever being looked at closely.)
    for (bool finished = false; !finished;) {
        int farthest = -1, seen_from = -1;
        double best = b.eps;
        for (int p = 0; p < count; ++p)
            for (size_t t = 0; t < b.triangles.size(); ++t)
                if (b.triangles[t].alive && b.height(b.triangles[t], p) > best) best = b.height(b.triangles[t], p), farthest = p, seen_from = static_cast<int>(t);
        if (farthest < 0) finished = true;
        else if (!b.absorb(farthest, seen_from)) return false;
    }

    std::vector<Triangle> surface;
    for (const Triangle& t : b.triangles)
        if (t.alive) surface.push_back(t);

    // Merge triangles that lie in one plane. Faces are grown one at a time from a seed triangle: a
    // neighbour joins if all three of its corners lie in the SEED's plane.
    // (Judging each triangle against its neighbour instead lets a face creep round a gentle curve, each
    // step within tolerance and the whole of it far from flat.)
    const int n = static_cast<int>(surface.size());
    auto neighbour_across = [&](int tri, int e) {  // the triangle on the other side of edge e of `tri`, or -1
        const int from = surface[static_cast<size_t>(tri)].v[e], to = surface[static_cast<size_t>(tri)].v[(e + 1) % 3];
        for (int other = 0; other < n; ++other)
            for (int k = 0; k < 3; ++k)
                if (surface[static_cast<size_t>(other)].v[k] == to && surface[static_cast<size_t>(other)].v[(k + 1) % 3] == from) return other;
        return -1;
    };
    std::vector<int> group(static_cast<size_t>(n), -1);  // the seed triangle of the face each triangle belongs to
    for (int seed = 0; seed < n; ++seed) {
        if (group[static_cast<size_t>(seed)] >= 0) continue;
        const Triangle& plane = surface[static_cast<size_t>(seed)];
        std::vector<int> frontier{seed};
        group[static_cast<size_t>(seed)] = seed;
        while (!frontier.empty()) {
            const int tri = frontier.back();
            frontier.pop_back();
            for (int e = 0; e < 3; ++e) {
                const int other = neighbour_across(tri, e);
                if (other < 0) return false;  // a hole in the surface
                if (group[static_cast<size_t>(other)] >= 0) continue;
                const Triangle& t = surface[static_cast<size_t>(other)];
                bool level = true;
                for (int k = 0; k < 3; ++k) level = level && std::fabs(b.height(plane, t.v[k])) <= b.eps;
                if (!level) continue;
                group[static_cast<size_t>(other)] = seed;
                frontier.push_back(other);
            }
        }
    }
    struct { const std::vector<int>& of; int find(int tri) const { return tri < 0 ? -1 : of[static_cast<size_t>(tri)]; } } groups{group};

    // Each group is one face. Its outline is made of the edges whose other side belongs to a different group;
    // follow them head to tail.
    std::vector<std::vector<int>> outlines;
    std::vector<int> seeds;
    for (int root = 0; root < n; ++root) {
        if (groups.find(root) != root) continue;
        std::vector<std::pair<int, int>> outline;
        for (int tri = 0; tri < n; ++tri) {
            if (groups.find(tri) != root) continue;
            for (int e = 0; e < 3; ++e)
                if (groups.find(neighbour_across(tri, e)) != root) outline.emplace_back(surface[static_cast<size_t>(tri)].v[e], surface[static_cast<size_t>(tri)].v[(e + 1) % 3]);
        }
        std::vector<int> loop{outline[0].first};
        for (int at = outline[0].second; at != outline[0].first;) {
            if (loop.size() > outline.size()) return false;  // not a single closed outline
            loop.push_back(at);
            const auto next = std::find_if(outline.begin(), outline.end(), [&](const std::pair<int, int>& e) { return e.first == at; });
            if (next == outline.end()) return false;
            at = next->second;
        }
        if (loop.size() != outline.size()) return false;
        outlines.push_back(loop);
        seeds.push_back(root);
    }

    // A real corner is where three or more faces meet. A point on only two is a point along the edge
    // between them (the triangles that made it a corner have all been merged away), so it goes, from both.
    for (bool removed = true; removed;) {
        removed = false;
        std::vector<int> faces_at(static_cast<size_t>(count), 0);
        for (const auto& loop : outlines)
            for (const int v : loop) ++faces_at[static_cast<size_t>(v)];
        for (auto& loop : outlines) removed = std::erase_if(loop, [&](int v) { return faces_at[static_cast<size_t>(v)] == 2; }) > 0 || removed;
    }

    std::vector<int> remap(static_cast<size_t>(count), -1);
    for (size_t k = 0; k < outlines.size(); ++k) {
        if (outlines[k].size() < 3) return false;
        // Merging by tolerance can, on a very thin cloud, produce an outline that is not convex. Refuse it.
        const std::vector<int>& loop = outlines[k];
        for (size_t i = 0; i < loop.size(); ++i) {
            const Vec from = b.at(loop[i]), side = b.at(loop[(i + 1) % loop.size()]) - from;
            Vec inward = cross(surface[static_cast<size_t>(seeds[k])].normal, side);
            inward = inward * (1 / inward.length());
            for (const int v : loop)
                if (dot(inward, b.at(v) - from) < -2 * b.eps) return false;
        }
        // The face's plane is its seed triangle's: a plane that by construction has every input point on or
        // behind it (to within the tolerance) and every corner of the face within the tolerance of it.
        Hull::Face face;
        const Triangle& plane = surface[static_cast<size_t>(seeds[k])];
        face.normal = {static_cast<Real>(plane.normal.x), static_cast<Real>(plane.normal.y), static_cast<Real>(plane.normal.z)};
        face.offset = static_cast<Real>(plane.offset);
        face.first = static_cast<int>(out.loops.size());
        face.count = static_cast<int>(outlines[k].size());
        for (const int v : outlines[k]) {
            if (remap[static_cast<size_t>(v)] < 0) {
                remap[static_cast<size_t>(v)] = static_cast<int>(out.vertices.size());
                out.vertices.push_back(points[v]);
            }
            out.loops.push_back(remap[static_cast<size_t>(v)]);
        }
        out.faces.push_back(face);
    }
    // Last checks that cost nothing. A solid has at least four faces, and a closed surface around one has
    // V - E + F = 2. Anything else means the cloud was too thin or too tangled for the tolerance, and no
    // hull is better than a wrong one.
    return out.faces.size() >= 4 && static_cast<int>(out.vertices.size()) - out.edge_count() + static_cast<int>(out.faces.size()) == 2;
}

}  // namespace

bool build_hull(const Vec3* points, int count, Hull& out) {
    out = Hull{};
    if (build(points, count, out)) return true;
    out = Hull{};
    return false;
}

}  // namespace phys3d
