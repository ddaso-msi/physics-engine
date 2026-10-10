#pragma once
// Stage 10, step 9: GJK and EPA, distance and penetration between ANY two convex shapes.
//
// The box and sphere routines in collision.cpp each know their shapes' geometry. These two algorithms know
// only one thing about a shape: its SUPPORT FUNCTION, "which of your points is farthest in direction d?".
// Anything convex can answer that, so one piece of code handles every pairing.
//
// The idea both rest on is the Minkowski difference  A - B = { a - b : a in A, b in B }.
//   - A and B overlap exactly when A - B contains the origin (some a equals some b).
//   - If they do not, the distance between them is the distance from the origin to A - B.
//   - If they do, the shortest way to pull them apart is the shortest way from the origin out of A - B.
// A - B is convex, and its support function is free: the farthest point of A - B along d is
// (farthest of A along d) - (farthest of B along -d). So we never build A - B; we only sample it.
//
// GJK (Gilbert, Johnson, Keerthi) finds the point of A - B closest to the origin. It keeps a SIMPLEX, up to
// four sampled points (a point, segment, triangle or tetrahedron), finds the point of the simplex closest to
// the origin, asks for the support point in the direction from there toward the origin, and repeats. It
// stops when the new support point is no nearer the origin than the simplex already gets (separated), or
// when the simplex encloses the origin (overlapping).
//
// EPA (the Expanding Polytope Algorithm) takes over in the second case. It starts from GJK's tetrahedron,
// which is inside A - B and around the origin, and grows it: find the face nearest the origin, ask for the
// support point along that face's normal, and if A - B reaches further that way, add the point and rebuild
// the faces around it. When the nearest face cannot be pushed out any more, it lies on the surface of A - B,
// and its normal and distance are the penetration direction and depth.
//
// One refinement, used here: a shape is a CORE plus a rounding radius. A sphere is a point with a radius, a
// capsule a segment with a radius. The algorithms run on the cores; if the cores are apart, the rounded
// shapes' distance is just the cores' distance minus the two radii, exact and with no EPA, even when the
// rounded shapes overlap. EPA runs only when the cores themselves overlap.
#include "body.hpp"
#include "hull.hpp"

namespace phys3d {

// A convex shape placed in the world, known only through its support function.
struct Convex {
    enum class Core {
        Point,    // with a radius: a sphere
        Segment,  // along the local y axis, from -half.y to +half.y; with a radius: a capsule
        Box,      // half extents `half`
        Points,   // the convex hull of a set of local points
    };

    Core core = Core::Point;
    Vec3 pos;
    Quat q;
    Vec3 half;
    const Vec3* points = nullptr;  // Core::Points: not owned; must outlive the Convex
    int count = 0;
    // Core::Points, optional: the same points as a built hull, which knows its faces. GJK and EPA do not
    // need it; contact patches do (a shape can only rest on a face it knows it has). Not owned either.
    const Hull* faces = nullptr;
    Real radius = 0;  // every point within this distance of the core belongs to the shape

    static Convex sphere(Real radius, Vec3 pos);
    static Convex capsule(Real half_length, Real radius, Vec3 pos, Quat q = {});
    static Convex box(Vec3 half, Vec3 pos, Quat q = {});
    static Convex hull(const Vec3* points, int count, Vec3 pos, Quat q = {}, Real radius = 0);
    // A built hull placed in the world: its corners for GJK and EPA, its faces for contact patches.
    static Convex hull(const Hull& built, Vec3 pos, Quat q = {}, Real radius = 0);
    static Convex of(const Body& body);

    // The core's point farthest along world direction d (any length, not zero).
    Vec3 core_support(Vec3 d) const;
    // The whole shape's.
    Vec3 support(Vec3 d) const { return core_support(d) + d.normalized() * radius; }
};

struct ClosestResult {
    // Surface to surface. Positive: apart by this much. Negative: overlapping, and -distance is the
    // penetration depth, the shortest move that separates them.
    Real distance = 0;
    // Unit vector from a toward b: the direction to move b (by -distance) to separate an overlapping pair,
    // or the direction of the gap between a separated one.
    Vec3 normal;
    // Apart: the closest point of each surface. Overlapping: the deepest point of each, so that
    // point_b + normal * depth == point_a.
    Vec3 point_a, point_b;
    int gjk_iterations = 0;
    int epa_iterations = 0;  // 0 unless the cores overlapped
    // False if an iteration cap or a degenerate polytope stopped the search early; the answer is then the
    // best found so far.
    bool converged = true;
};

ClosestResult closest(const Convex& a, const Convex& b);

// The two small geometric routines GJK is built from, exposed so they can be tested on their own.
namespace detail {
// The point of segment ab closest to the ORIGIN, as weights: the point is a * out[0] + b * out[1].
void closest_on_segment(Vec3 a, Vec3 b, Real out[2]);
// The same for triangle abc.
void closest_on_triangle(Vec3 a, Vec3 b, Vec3 c, Real out[3]);
}  // namespace detail

}  // namespace phys3d
