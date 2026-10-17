#pragma once
// Stage 10, step 23: joints solved exactly, a tree of them at a time. Internal to the engine.
//
// Gauss-Seidel passes a correction along a chain one joint per sweep, so ten sweeps leave a chain of ten welded
// links soft: a welded rod 7 m long sagged and opened by 2.3 mm on average and 14 mm at worst. For joints that
// form a tree (chains, ragdolls, anything without a loop) the equations can instead be solved exactly, in time
// proportional to their number (Baraff, "Linear-Time Dynamics using Lagrange Multipliers", 1996). Take the
// bodies' velocity changes dv and the rows' impulses as unknowns together:
//
//     [ M  J^T ] [ dv ]   [ 0 ]            M   the bodies' masses and inertias (6 x 6 each)
//     [ J   0  ] [ mu ] = [ c ]            J   the rows, each touching one or two bodies
//                                          c   -(J v + bias): what each row's speed must change by
//
// The first line says dv = -M^-1 J^T mu (so the impulse is -mu), the second that the rows end up satisfied. The
// matrix has a block for each body and each joint, and off the diagonal a block for each joint and body it
// holds: drawn as a graph, the bodies and joints with an edge for each block, it is the joints' own tree. Gaussian
// elimination taken leaves first then creates no new blocks: eliminating a leaf changes only its parent's
// diagonal block. So the factoring costs a few 6 x 6 products per body and joint, and so does each solve.
//
// Only the rows that hold both ways, rigidly (anchors, hinge and weld axes, slider rails) are solved here. Limits,
// motors, the cone and springs are inequalities or soft, and are left to Gauss-Seidel; so are joints that are
// not part of a tree (a loop, or a structure fixed to the world in two places, which is a loop through the
// world), and joints of a body that cannot turn about some axis (its mass matrix has no inverse). Gauss-Seidel
// sweeps every row, in the trees too; the trees' exact solve comes at the end of the last sweep.
#include <phys3d/joints.hpp>

#include <cstdint>
#include <vector>

namespace phys3d {

class JointTrees {
public:
    // Finds the trees among `rows` (built for the velocity pass; a joint's rows are rows[first[k] ..
    // first[k + 1]), the last entry being rows.size()) and factors them for the bodies' present poses.
    void build(const std::vector<Body>& bodies, const std::vector<JointRow>& rows, const std::vector<int>& first);
    // One exact solve of the tree rows from the bodies' current velocities: changes the velocities and adds to
    // the rows' accumulated impulses.
    void solve(std::vector<JointRow>& rows);

private:
    // (The arrays are left uninitialised, by the empty constructors: each node is a kilobyte, and clearing two
    // thousand of them every step cost more than factoring them.)
    struct Node {
        Node() {}
        int size;          // 6 for a body, its number of rows for a joint
        int parent;        // node index; -1 for a tree's root
        Body* body;        // a body node's body; null for a joint
        int rows[6];       // a joint node's rows
        double Dinv[6][6];  // the inverse of this node's diagonal block once its children are eliminated
        double T[6][6];     // Dinv times the block coupling this node to its parent (size x parent's size)
    };
    // What a solve needs of a node, packed: its blocks in single precision, only as large as they are, one node
    // after another in elimination order. (The factoring is done in double precision.)
    struct Step {
        int size, parent, parent_size;  // parent: index of the parent's Step, -1 for a root
        int data;                       // where Dinv (size x size) and then T (size x parent_size) start in data_
        int vec;                        // where y and x (size each) start in vec_
        Body* body;                     // null for a joint
        int rows[6];
    };
    // Only while building: a joint node's body nodes, and its rows as they act on each (row, then linear and
    // angular), from which the blocks coupling it to its neighbours are taken.
    struct Links {
        Links() {}
        int node[2];  // -1: held to the world or a fixed body
        double J[2][6][6];
    };
    int add_node(int size, Body* body);  // returns its index
    // The block coupling child c to its parent: H[c][parent], size(c) x size(parent).
    void child_block(int c, double out[6][6]) const;

    std::vector<Node> nodes_;   // the first used_ are this step's
    std::vector<Links> links_;  // per node (unused for a body)
    int used_ = 0;
    std::vector<Step> steps_, sorted_;  // every tree, packed for solving, and ordered by height (see build)
    std::vector<float> data_, vec_;
    std::vector<int> body_node_, adjacency_, adjacent_start_, piece_, post_, step_of_;
    std::vector<std::uint8_t> seen_;
    std::vector<std::pair<int, std::size_t>> walk_;
};

}  // namespace phys3d
