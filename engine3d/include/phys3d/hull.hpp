#pragma once
// Stage 10, step 12: building a convex hull.
//
// GJK and EPA are content with a convex shape given as a cloud of points: the farthest point in a direction
// is all they ask. A contact PATCH needs more. To rest on a face a shape has to know it has one: which
// points form it, in what order, facing which way. This turns the cloud into that: the corner points that
// matter (the rest are inside and dropped) and the flat polygons between them.
//
// The method is the one EPA uses to grow its polytope, run to completion (it is Quickhull without the
// bookkeeping that makes Quickhull fast on large inputs):
//   1. Start with a tetrahedron of four of the points, spread as widely as possible.
//   2. Find the point farthest outside the current surface. If none is outside, stop.
//   3. Remove every triangle that point can see. The hole left has a rim, the horizon. Close it with a fan of
//      new triangles from the rim to the point. Go to 2.
// That gives a surface of triangles. A cube's face comes out as two of them, so as a last step triangles
// that lie in one plane are merged into a single polygon, and the result is what the shape "is": six
// squares, not twelve triangles.
#include "math.hpp"

#include <vector>

namespace phys3d {

struct Hull {
    struct Face {
        Vec3 normal;      // unit, pointing out of the shape
        Real offset = 0;  // the face lies in the plane dot(normal, p) == offset
        int first = 0;    // its corners are loops[first .. first + count), counter-clockwise seen from outside
        int count = 0;
    };

    std::vector<Vec3> vertices;  // the corners of the hull, and nothing else
    std::vector<Face> faces;
    std::vector<int> loops;      // indices into `vertices`, face after face

    int corner(const Face& f, int i) const { return loops[static_cast<size_t>(f.first + ((i % f.count) + f.count) % f.count)]; }
    const Vec3& point(const Face& f, int i) const { return vertices[static_cast<size_t>(corner(f, i))]; }
    // Every edge is shared by two faces, so there are half as many edges as there are face corners.
    int edge_count() const { return static_cast<int>(loops.size()) / 2; }
};

// Builds the hull of `count` points. Points inside the hull, on its faces or along its edges do not become
// vertices. Returns false, leaving `out` empty, if the points do not enclose a volume (fewer than four, or
// all in one plane, to within a hundred-thousandth of their spread).
bool build_hull(const Vec3* points, int count, Hull& out);

}  // namespace phys3d
