#include "test.hpp"
#include <phys3d/gjk.hpp>
#include <phys3d/hull.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

using namespace phys3d;

namespace {

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Vec3 in_box(Real r) { return {range(-r, r), range(-r, r), range(-r, r)}; }
    Vec3 direction() {
        for (;;) {
            const Vec3 v = in_box(1);
            if (v.length_sq() > 0.01f && v.length_sq() <= 1) return v.normalized();
        }
    }
    Quat rotation() { return Quat::from_axis_angle(direction(), range(-kPi, kPi)); }
};

Real spread(const std::vector<Vec3>& points) {
    Vec3 lo = points[0], hi = points[0];
    for (const Vec3& p : points)
        for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], p[k]), hi[k] = std::max(hi[k], p[k]);
    return distance(lo, hi);
}

// Checks everything a hull must be, against the points it was built from. Returns what is wrong, or nullptr.
const char* audit(const Hull& h, const std::vector<Vec3>& points) {
    const Real tol = 3e-5f * spread(points);
    if (h.faces.size() < 4 || h.vertices.size() < 4) return "fewer than four faces or vertices";
    std::set<std::pair<int, int>> edges;
    std::vector<int> uses(h.vertices.size(), 0);
    for (const Hull::Face& f : h.faces) {
        if (f.count < 3) return "a face with fewer than three corners";
        Real turning = 0;  // positive if the corners go round counter-clockwise
        if (std::fabs(f.normal.length() - 1) > 1e-4f) return "a normal that is not a unit vector";
        for (int i = 0; i < f.count; ++i) {
            // Flat: every corner in the face's plane.
            if (std::fabs(dot(f.normal, h.point(f, i)) - f.offset) > tol) return "a face that is not flat";
            // Convex and counter-clockwise seen from outside: every other corner is to the left of every
            // side (give or take the tolerance: a face merged from triangles that were only nearly level
            // can have a side that is out of line by that much).
            const Vec3 side = h.point(f, i + 1) - h.point(f, i), inward = cross(f.normal, side).normalized();
            for (int j = 0; j < f.count; ++j)
                if (dot(inward, h.point(f, j) - h.point(f, i)) < -tol) return "a face that is not convex and counter-clockwise";
            turning += dot(cross(side, h.point(f, i + 2) - h.point(f, i + 1)), f.normal);
            if (!edges.insert({h.corner(f, i), h.corner(f, i + 1)}).second) return "an edge used twice in the same direction";
            ++uses[static_cast<size_t>(h.corner(f, i))];
        }
        if (turning <= 0) return "a face wound the wrong way";
        // A supporting plane: nothing of the input is in front of it.
        for (const Vec3& p : points)
            if (dot(f.normal, p) - f.offset > tol) return "an input point outside a face";
    }
    // Closed: every edge has its twin, running the other way, on the neighbouring face.
    for (const auto& e : edges)
        if (!edges.count({e.second, e.first})) return "an edge with no face on its other side";
    // Euler's formula for anything shaped like a ball: V - E + F = 2.
    if (static_cast<int>(h.vertices.size()) - h.edge_count() + static_cast<int>(h.faces.size()) != 2) return "V - E + F is not 2";
    for (size_t v = 0; v < h.vertices.size(); ++v) {
        if (uses[v] < 3) return "a vertex on fewer than three faces";
        if (std::find_if(points.begin(), points.end(), [&](const Vec3& p) { return distance(p, h.vertices[v]) == 0; }) == points.end())
            return "a vertex that is not one of the input points";
    }
    // Merged: no two neighbouring faces lie squarely in the same plane (every corner of each within a third
    // of the builder's own tolerance of the other's plane; closer calls than that could go either way).
    const Real level = 3e-6f * spread(points);
    for (const Hull::Face& f : h.faces)
        for (const Hull::Face& g : h.faces) {
            if (&f == &g || dot(f.normal, g.normal) <= 0) continue;
            bool neighbours = false, coplanar = true;
            for (int i = 0; i < f.count; ++i) {
                coplanar = coplanar && std::fabs(dot(g.normal, h.point(f, i)) - g.offset) <= level;
                for (int j = 0; j < g.count; ++j) neighbours = neighbours || (h.corner(f, i) == h.corner(g, j + 1) && h.corner(f, i + 1) == h.corner(g, j));
            }
            for (int j = 0; j < g.count; ++j) coplanar = coplanar && std::fabs(dot(f.normal, h.point(g, j)) - f.offset) <= level;
            if (neighbours && coplanar) return "two coplanar faces left unmerged";
        }
    return nullptr;
}

Hull hull_of(const std::vector<Vec3>& points) {
    Hull h;
    CHECK(build_hull(points.data(), static_cast<int>(points.size()), h));
    const char* wrong = audit(h, points);
    if (wrong) std::printf("    %s (%zu points -> %zu vertices, %zu faces)\n", wrong, points.size(), h.vertices.size(), h.faces.size());
    CHECK(wrong == nullptr);
    return h;
}

std::vector<Vec3> cube_corners(Vec3 half = {1, 1, 1}) {
    std::vector<Vec3> p;
    for (int i = 0; i < 8; ++i) p.push_back({(i & 1) ? half.x : -half.x, (i & 2) ? half.y : -half.y, (i & 4) ? half.z : -half.z});
    return p;
}
// How many faces have `corners` corners.
int faces_with(const Hull& h, int corners) {
    return static_cast<int>(std::count_if(h.faces.begin(), h.faces.end(), [&](const Hull::Face& f) { return f.count == corners; }));
}
// The volume enclosed, by the divergence theorem: each face, fanned into triangles, makes tetrahedra with
// the origin, and their signed volumes add up to the volume of the solid.
Real volume(const Hull& h) {
    Real six = 0;
    for (const Hull::Face& f : h.faces)
        for (int i = 1; i + 1 < f.count; ++i) six += dot(h.point(f, 0), cross(h.point(f, i), h.point(f, i + 1)));
    return six / 6;
}

}  // namespace

// ---- shapes whose hull is known -------------------------------------------------------------------

TEST(the_hull_of_a_cubes_corners_is_six_squares) {
    const Hull h = hull_of(cube_corners({1, 2, 3}));
    CHECK(h.vertices.size() == 8);
    CHECK(h.faces.size() == 6 && faces_with(h, 4) == 6);
    CHECK(h.edge_count() == 12);
    CHECK_NEAR(volume(h), 48, 1e-3);
    // One face for each of the six directions, at the right distance.
    for (const Vec3& n : {Vec3{1, 0, 0}, Vec3{-1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, -1, 0}, Vec3{0, 0, 1}, Vec3{0, 0, -1}}) {
        const auto face = std::find_if(h.faces.begin(), h.faces.end(), [&](const Hull::Face& f) { return distance(f.normal, n) < 1e-5f; });
        CHECK(face != h.faces.end());
        if (face != h.faces.end()) CHECK_NEAR(face->offset, std::fabs(n.x) * 1 + std::fabs(n.y) * 2 + std::fabs(n.z) * 3, 1e-5);
    }
}

TEST(points_inside_on_faces_and_along_edges_do_not_become_vertices) {
    std::vector<Vec3> points = cube_corners();
    Lcg rng;
    for (int i = 0; i < 60; ++i) points.push_back(rng.in_box(0.99f));                       // inside
    for (int i = 0; i < 20; ++i) points.push_back({1, rng.range(-0.9f, 0.9f), rng.range(-0.9f, 0.9f)});  // on a face
    for (int i = 0; i < 10; ++i) points.push_back({1, 1, rng.range(-0.9f, 0.9f)});          // along an edge
    points.push_back({0, 0, 1});    // the middle of a face
    points.push_back({1, -1, 0});   // the middle of an edge
    for (int i = 0; i < 8; ++i) points.push_back(points[static_cast<size_t>(i)]);  // every corner twice
    // Whatever order they come in.
    for (int shuffle = 0; shuffle < 30; ++shuffle) {
        for (size_t i = points.size() - 1; i > 0; --i) std::swap(points[i], points[static_cast<size_t>(rng.next() * static_cast<Real>(i + 1))]);
        const Hull h = hull_of(points);
        CHECK(h.vertices.size() == 8);
        CHECK(h.faces.size() == 6 && faces_with(h, 4) == 6);
    }
}

TEST(a_lattice_of_points_has_the_hull_of_its_outer_box) {
    // 4 x 3 x 5 points on a grid: nearly every one is level with some face or in line with some edge.
    std::vector<Vec3> points;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 5; ++k) points.push_back({static_cast<Real>(i) * 0.5f, static_cast<Real>(j) * 0.7f, static_cast<Real>(k) * 0.3f});
    Lcg rng;
    for (int trial = 0; trial < 40; ++trial) {
        const Quat q = trial ? rng.rotation() : Quat{};
        const Vec3 shift = trial ? rng.in_box(2) : Vec3{};
        std::vector<Vec3> turned;
        for (const Vec3& p : points) turned.push_back(rotate(q, p) + shift);
        const Hull h = hull_of(turned);
        CHECK(h.vertices.size() == 8);
        CHECK(h.faces.size() == 6 && faces_with(h, 4) == 6);
        CHECK_NEAR(volume(h), 1.5 * 1.4 * 1.2, 2e-3);
    }
}

TEST(the_hulls_of_the_regular_solids_have_the_right_faces) {
    const Hull tetra = hull_of({{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}});
    CHECK(tetra.vertices.size() == 4 && tetra.faces.size() == 4 && faces_with(tetra, 3) == 4);

    const Hull octa = hull_of({{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}});
    CHECK(octa.vertices.size() == 6 && octa.faces.size() == 8 && faces_with(octa, 3) == 8);
    CHECK_NEAR(volume(octa), 4.0 / 3, 1e-4);

    // Icosahedron: the twelve points (0, +-1, +-g), (+-1, +-g, 0), (+-g, 0, +-1), g the golden ratio.
    const Real g = (1 + std::sqrt(5.0f)) / 2;
    std::vector<Vec3> twelve;
    for (const Real a : {-1.0f, 1.0f})
        for (const Real b : {-g, g}) {
            twelve.push_back({0, a, b});
            twelve.push_back({a, b, 0});
            twelve.push_back({b, 0, a});
        }
    const Hull icosa = hull_of(twelve);
    CHECK(icosa.vertices.size() == 12 && icosa.faces.size() == 20 && faces_with(icosa, 3) == 20);
    CHECK(icosa.edge_count() == 30);
}

TEST(a_prism_keeps_its_many_sided_ends_as_single_faces) {
    for (int sides : {3, 5, 6, 8, 12}) {
        std::vector<Vec3> points;
        for (int i = 0; i < sides; ++i) {
            const Real a = 2 * kPi * static_cast<Real>(i) / static_cast<Real>(sides);
            points.push_back({std::cos(a), -0.5f, std::sin(a)});
            points.push_back({std::cos(a), 0.5f, std::sin(a)});
        }
        const Hull h = hull_of(points);
        CHECK(static_cast<int>(h.vertices.size()) == 2 * sides);
        CHECK(static_cast<int>(h.faces.size()) == sides + 2);
        CHECK(faces_with(h, sides) == (sides == 4 ? 6 : 2));
        CHECK(faces_with(h, 4) == sides);
        // The area of a regular polygon inscribed in the unit circle, times the height 1.
        CHECK_NEAR(volume(h), 0.5 * sides * std::sin(2 * kPi / static_cast<Real>(sides)), 1e-4);
    }
}

// ---- any cloud of points -------------------------------------------------------------------------

TEST(points_on_a_sphere_are_all_vertices_and_all_faces_are_triangles) {
    Lcg rng;
    for (int n : {4, 5, 8, 20, 60, 150}) {
        std::vector<Vec3> points;
        for (int i = 0; i < n; ++i) points.push_back(rng.direction() * 2.0f);
        const Hull h = hull_of(points);
        CHECK(static_cast<int>(h.vertices.size()) == n);
        // If every face is a triangle, Euler's formula fixes the counts: F = 2V - 4, E = 3V - 6. Each pair of
        // triangles merged into a quadrilateral takes one off both.
        const int quads = faces_with(h, 4);
        CHECK(static_cast<int>(h.faces.size()) == 2 * n - 4 - quads && faces_with(h, 3) == 2 * n - 4 - 2 * quads);
        CHECK(h.edge_count() == 3 * n - 6 - quads);
        // Four random points on a sphere all but never share a plane; with a tolerance, a few among
        // hundreds of neighbouring pairs can come close enough.
        CHECK(quads <= n / 50);
    }
}

TEST(random_clouds_give_valid_hulls_with_the_same_support_function) {
    Lcg rng;
    for (int trial = 0; trial < 400; ++trial) {
        const int n = 4 + static_cast<int>(rng.next() * 60);
        const Real scale = trial % 3 == 0 ? 0.01f : trial % 3 == 1 ? 1.0f : 100.0f;
        const Vec3 shift = trial % 5 == 0 ? rng.in_box(50) * scale : Vec3{};
        std::vector<Vec3> points;
        for (int i = 0; i < n; ++i) points.push_back(rng.in_box(1) * scale + shift);
        const Hull h = hull_of(points);
        // The hull reaches exactly as far as the cloud in every direction: it has lost nothing that mattered.
        for (int k = 0; k < 60; ++k) {
            const Vec3 d = rng.direction();
            Real cloud = -1e30f, corners = -1e30f;
            for (const Vec3& p : points) cloud = std::max(cloud, dot(p, d));
            for (const Vec3& p : h.vertices) corners = std::max(corners, dot(p, d));
            CHECK(cloud == corners);
        }
    }
}

TEST(the_faces_enclose_exactly_what_gjk_says_the_cloud_encloses) {
    // An independent opinion on "inside": GJK on the raw points, which never sees the faces.
    Lcg rng;
    int inside = 0, outside = 0;
    for (int trial = 0; trial < 60; ++trial) {
        std::vector<Vec3> points;
        const int n = 5 + static_cast<int>(rng.next() * 30);
        for (int i = 0; i < n; ++i) points.push_back(rng.in_box(1));
        const Hull h = hull_of(points);
        const Convex cloud = Convex::hull(points.data(), n, {0, 0, 0});
        for (int k = 0; k < 60; ++k) {
            const Vec3 q = rng.in_box(1.1f);
            Real deepest_plane = -1e30f;  // how far q is in front of the face it is most in front of
            for (const Hull::Face& f : h.faces) deepest_plane = std::max(deepest_plane, dot(f.normal, q) - f.offset);
            const Real by_gjk = closest(cloud, Convex::sphere(0, q)).distance;
            if (std::fabs(deepest_plane) < 2e-3f) continue;  // on the surface: too close to call
            CHECK((deepest_plane < 0) == (by_gjk < 0));
            (deepest_plane < 0 ? inside : outside) += 1;
        }
    }
    CHECK(inside > 500 && outside > 500);
}

TEST(a_hull_does_not_depend_on_the_order_of_the_points) {
    Lcg rng;
    for (int trial = 0; trial < 50; ++trial) {
        std::vector<Vec3> points;
        for (int i = 0; i < 30; ++i) points.push_back(rng.in_box(1));
        const Hull first = hull_of(points);
        for (size_t i = points.size() - 1; i > 0; --i) std::swap(points[i], points[static_cast<size_t>(rng.next() * static_cast<Real>(i + 1))]);
        const Hull second = hull_of(points);
        CHECK(first.vertices.size() == second.vertices.size() && first.faces.size() == second.faces.size());
        CHECK_NEAR(volume(first), volume(second), 1e-4);
    }
}

// ---- thin shapes, where the tolerance is tested -----------------------------------------------------

namespace {

// A cloud a thousand times wider than it is thick: every triangle of its top and bottom is close to level
// with its neighbours, and the builder's tolerance is a hundredth of the thickness.
std::vector<Vec3> flat_cloud(int trial) {
    Lcg rng;
    rng.s = static_cast<std::uint32_t>(trial) * 2654435761u + 99;
    const int n = 8 + static_cast<int>(rng.next() * 55);
    const Real scale = trial % 3 == 0 ? 0.01f : trial % 3 == 1 ? 1.0f : 100.0f;
    std::vector<Vec3> points;
    for (int i = 0; i < n; ++i) {
        Vec3 p = rng.in_box(1);
        p.z *= 1e-3f;
        points.push_back(p * scale);
    }
    return points;
}

}  // namespace

TEST(very_flat_clouds_get_a_valid_hull_or_none) {
    int built = 0;
    for (int trial = 0; trial < 600; ++trial) {
        const std::vector<Vec3> points = flat_cloud(trial);
        Hull h;
        if (!build_hull(points.data(), static_cast<int>(points.size()), h)) continue;
        ++built;
        const char* wrong = audit(h, points);
        if (wrong) std::printf("    flat cloud %d: %s\n", trial, wrong);
        CHECK(wrong == nullptr);
    }
    CHECK(built > 590);  // refusing is allowed, but it must be rare
}

TEST(flat_clouds_whose_merged_faces_would_not_be_convex_are_refused_or_built_right) {
    // Found by a stress run: on these the tolerance lets triangles join a face whose outline then has a dent.
    int refused = 0;
    for (int trial : {1944, 3952, 6566, 6896, 10729, 21895, 23573, 25028, 33420, 35304, 40204, 41532}) {
        const std::vector<Vec3> points = flat_cloud(trial);
        Hull h;
        if (!build_hull(points.data(), static_cast<int>(points.size()), h)) {
            ++refused;
            CHECK(h.vertices.empty() && h.faces.empty() && h.loops.empty());  // refused late, and still left empty
            continue;
        }
        CHECK(audit(h, points) == nullptr);
    }
    CHECK(refused > 0);
}

// ---- clouds with no hull -------------------------------------------------------------------------

TEST(points_that_enclose_no_volume_are_refused) {
    Hull h;
    h.vertices.push_back({1, 2, 3});  // whatever was there before is cleared
    const Vec3 three[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    CHECK(!build_hull(three, 3, h));
    CHECK(h.vertices.empty() && h.faces.empty() && h.loops.empty());

    Lcg rng;
    std::vector<Vec3> flat, line, same;
    const Quat q = rng.rotation();
    for (int i = 0; i < 40; ++i) {
        flat.push_back(rotate(q, {rng.range(-1, 1), rng.range(-1, 1), 0.3f}));
        line.push_back(rotate(q, {rng.range(-1, 1), 0.2f, 0.3f}));
        same.push_back({0.5f, -0.25f, 2});
    }
    CHECK(!build_hull(flat.data(), 40, h));
    CHECK(!build_hull(line.data(), 40, h));
    CHECK(!build_hull(same.data(), 40, h));
    CHECK(h.vertices.empty() && h.faces.empty());

    // A cloud with every point at the same x is still a perfectly good solid.
    std::vector<Vec3> slab;
    for (int i = 0; i < 4; ++i) slab.push_back({2, (i & 1) ? 1.0f : -1.0f, (i & 2) ? 1.0f : -1.0f});
    slab.push_back({2, 0, 0});
    CHECK(!build_hull(slab.data(), 5, h));  // ...unless they really are all in one plane
    slab.push_back({3, 0, 0});
    CHECK(build_hull(slab.data(), 6, h));
    CHECK(h.vertices.size() == 5 && h.faces.size() == 5);  // a pyramid on a square
}
