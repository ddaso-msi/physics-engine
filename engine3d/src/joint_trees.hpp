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
// motors, the cone and springs are inequalities or soft, and are left to Gauss-Seidel; so are joints of a body that
// cannot turn about some axis (its mass matrix has no inverse). Gauss-Seidel sweeps every row, in the trees too; the
// trees' exact solve comes at the end of the last sweep.
//
// Except in the position pass, where a limit or cone that is past its limit is taken in with the rest (step 31): it
// has an error to remove, as they do, and removing them one after the other undid each other. Turning a thigh back
// inside its cone swings the hip's anchor open, and closing the anchor again (the tree's answer, the least motion
// that closes it) turns the thigh back out: a ragdoll's hip went round that circle four times a step, 6% closer each
// time. Solved together, both are met at once. (A joint's rows still number at most six: a ball joint's three and its
// cone; a hinge's five and one of its two limits, which cannot both be passed at once; a slider's likewise.)
//
// Loops (step 27): a ring of joints, or a structure fixed to the world in two places (a loop through the world), is
// not a tree. Each loop is cut at one joint, chosen as the joints are taken in turn: a joint that would join two
// bodies already joined (or a body to the world when it already is) is a cut joint, and the rest are trees. The cut
// joints' rows are then put back exactly, through the tree: an impulse lambda on them moves the tree's bodies by
// A lambda, A found by one tree solve per cut row with that row's impulse on the right (so the tree's own rows still
// hold), and the cut rows' speeds by W lambda, W = J_cut A. After each tree solve, lambda = W^-1 (what the cut rows
// still need), and A lambda is applied. A tree with more than kMaxCutRows of them is left to Gauss-Seidel.
#include <phys3d/joints.hpp>

#include <cstdint>
#include <vector>

namespace phys3d {

class JointTrees {
public:
    // The numbers the solves work with. Double: in single precision the answer's body part and its impulses part
    // disagreed by a few parts in ten thousand over a chain of light links, and either way something was off: the
    // velocities with the impulses kept, or the rows with exactness.
    using Scalar = double;

    // Finds the trees among `rows` (built for the velocity pass; a joint's rows are rows[first[k] ..
    // first[k + 1]), the last entry being rows.size()) and factors them for the bodies' present poses.
    // If `taken` is given it receives, for each joint k, whether its rows were taken into a tree. With `violated`
    // (for the position pass) a one-sided row that is past its limit is taken too, as if it held both ways, and
    // its joint's entry in `taken` is 2 rather than 1.
    void build(const std::vector<Body>& bodies, const std::vector<JointRow>& rows, const std::vector<int>& first,
               std::vector<std::uint8_t>* taken = nullptr, bool violated = false);
    // One exact solve of the tree rows from the bodies' current velocities: changes the velocities and adds to
    // the rows' accumulated impulses.
    void solve(std::vector<JointRow>& rows) { pass(rows, false, 0, 0, nullptr); }
    // The same solve for the rows' position errors (rows built for the position pass), each limited to the
    // most a row may correct at once: moves and turns the bodies so as to remove them, to first order.
    // `largest[k]` receives, for each joint k taken into a tree, the largest error any row of its tree had
    // (metres or radians).
    void correct(std::vector<JointRow>& rows, Real max_linear, Real max_angular, std::vector<Real>& largest) {
        pass(rows, true, max_linear, max_angular, &largest);
    }

private:
    // (The arrays are left uninitialised, by the empty constructors: each node is a kilobyte, and clearing two
    // thousand of them every step cost more than factoring them.)
    struct Node {
        Node() {}
        int size;          // 6 for a body, its number of rows for a joint
        int parent;        // node index; -1 for a tree's root
        Body* body;        // a body node's body; null for a joint
        int rows[6];       // a joint node's rows
        int joint;         // a joint node's joint (its index k in `first`)
        double Dinv[6][6];  // the inverse of this node's diagonal block once its children are eliminated
        double T[6][6];     // Dinv times the block coupling this node to its parent (size x parent's size)
    };
    // What a solve needs of a node, packed: its blocks, only as large as they are, one node after another in
    // elimination order.
    struct Step {
        int size, parent, parent_size;  // parent: index of the parent's Step, -1 for a root
        int data;                       // where Dinv (size x size) and then T (size x parent_size) start in data_
        int vec;                        // where y and x (size each) start in vec_
        Body* body;                     // null for a joint
        int rows[6];
        int tree;                       // which tree it is in
        int joint;                      // a joint's index k in `first`; -1 for a body
    };
    // Only while building: a joint node's body nodes, and its rows as they act on each (row, then linear and
    // angular), from which the blocks coupling it to its neighbours are taken.
    struct Links {
        Links() {}
        int node[2];  // -1: held to the world or a fixed body
        double J[2][6][6];
    };
    // A joint that would close a loop: its two-way rows, and its bodies (index into the bodies; -1 for fixed).
    struct Cut {
        int joint, size, rows[6], body[2];
    };
    // A tree's cut rows, put back through it (see above).
    static constexpr int kMaxCutRows = 24;
    struct Loop {
        int tree = 0, size = 0;
        int rows[kMaxCutRows];
        int joint[kMaxCutRows];    // the joint (index k in `first`) each row belongs to
        int step[kMaxCutRows][2];  // the Step of each row's bodies, -1 for a fixed one
        Scalar J[kMaxCutRows][2][6];
        size_t response = 0;  // where A starts in response_: column by column, each the tree's x in its step order
        size_t length = 0;    // floats in one column
        std::vector<double> W;  // the cut rows' speeds answer an impulse lambda on them by W lambda
        std::vector<double> L;  // N = L L^T, N = (W W) over the rows not dropped (see close_loops)
        std::uint8_t dropped[kMaxCutRows];  // rows the others (with the tree) already imply: given no impulse
    };
    void solve_tree(int tree);  // one solve of a single tree, from the y already in vec_, leaving x there
    void close_loops(std::vector<JointRow>& rows, bool positions, Real max_linear, Real max_angular);
    int add_node(int size, Body* body);  // returns its index
    void pass(std::vector<JointRow>& rows, bool positions, Real max_linear, Real max_angular, std::vector<Real>* largest);
    // The block coupling child c to its parent: H[c][parent], size(c) x size(parent).
    void child_block(int c, double out[6][6]) const;

    std::vector<Node> nodes_;   // the first used_ are this step's
    std::vector<Links> links_;  // per node (unused for a body)
    int used_ = 0;
    std::vector<Step> steps_, sorted_;  // every tree, packed for solving, and ordered by height (see build)
    int trees_ = 0;
    std::vector<Real> tree_error_;
    std::vector<int> union_, tree_start_, tree_steps_, body_step_;  // tree_steps_[tree_start_[t] ..]: tree t's steps, in order
    std::vector<Cut> cuts_;
    std::vector<Loop> loops_;
    std::vector<Scalar> response_;
    std::vector<Scalar> data_, vec_;
    std::vector<int> body_node_, adjacency_, adjacent_start_, piece_, post_, step_of_;
    std::vector<std::uint8_t> seen_;
    std::vector<std::pair<int, std::size_t>> walk_;
};

}  // namespace phys3d
