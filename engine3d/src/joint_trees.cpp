#include "joint_trees.hpp"

#include <algorithm>
#include <array>
#include <utility>
#include <cmath>
#include <limits>
#include <span>

namespace phys3d {

namespace {

constexpr Real kInf = std::numeric_limits<Real>::infinity();

// Inverts in place the symmetric n x n matrix m, which is positive definite (sign 1: a body's block) or negative
// definite (sign -1: a joint's). By Cholesky: sign m = L L^T, so m^-1 = sign L^-T L^-1. False if it is not
// definite, to a billionth of its largest diagonal entry: rows that are not independent. (n is a template
// argument so the loops are unrolled; the factoring is otherwise most of the trees' cost.)
template <int n>
bool invert(double m[6][6], double sign) {
    double L[6][6], scale = 0;
    for (int i = 0; i < n; ++i) scale = std::max(scale, sign * m[i][i]);
    for (int c = 0; c < n; ++c) {
        double d = sign * m[c][c];
        for (int k = 0; k < c; ++k) d -= L[c][k] * L[c][k];
        if (!(d > 1e-9 * scale)) return false;
        L[c][c] = std::sqrt(d);
        for (int i = c + 1; i < n; ++i) {
            double e = sign * m[i][c];
            for (int k = 0; k < c; ++k) e -= L[i][k] * L[c][k];
            L[i][c] = e / L[c][c];
        }
    }
    double Li[6][6];  // L^-1, lower triangular
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < j; ++i) Li[i][j] = 0;
        Li[j][j] = 1 / L[j][j];
        for (int i = j + 1; i < n; ++i) {
            double sum = 0;
            for (int k = j; k < i; ++k) sum -= L[i][k] * Li[k][j];
            Li[i][j] = sum / L[i][i];
        }
    }
    for (int i = 0; i < n; ++i)  // m^-1 = sign Li^T Li, symmetric
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = i; k < n; ++k) sum += Li[k][i] * Li[k][j];
            m[i][j] = m[j][i] = sign * sum;
        }
    return true;
}

// Eliminating a child (size N) into its parent (size P): T = Dinv B, kept for the solves, and the parent's block
// less B^T Dinv B = B^T T (symmetric). B is the block coupling child to parent, N x P.
template <int N, int P>
void eliminate(const double Dinv[6][6], const double B[6][6], double T[6][6], double D[6][6]) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < P; ++j) {
            double sum = 0;
            for (int k = 0; k < N; ++k) sum += Dinv[i][k] * B[k][j];
            T[i][j] = sum;
        }
    for (int i = 0; i < P; ++i)
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = 0; k < N; ++k) sum += B[k][i] * T[k][j];
            D[i][j] -= sum;
            if (j != i) D[j][i] -= sum;
        }
}
using Invert = bool (*)(double[6][6], double);
using Eliminate = void (*)(const double[6][6], const double[6][6], double[6][6], double[6][6]);
constexpr std::array<Invert, 7> kInvert = {invert<1>, invert<1>, invert<2>, invert<3>, invert<4>, invert<5>, invert<6>};
template <int N, int... P>
constexpr std::array<Eliminate, 7> eliminations(std::integer_sequence<int, P...>) { return {eliminate<N, P>...}; }
constexpr std::array<std::array<Eliminate, 7>, 7> kEliminate = {
    eliminations<1>(std::make_integer_sequence<int, 7>{}), eliminations<1>(std::make_integer_sequence<int, 7>{}),
    eliminations<2>(std::make_integer_sequence<int, 7>{}), eliminations<3>(std::make_integer_sequence<int, 7>{}),
    eliminations<4>(std::make_integer_sequence<int, 7>{}), eliminations<5>(std::make_integer_sequence<int, 7>{}),
    eliminations<6>(std::make_integer_sequence<int, 7>{})};

// Three of the solve's numbers as a vector of the engine's.
Vec3 as_vec3(const JointTrees::Scalar* p) { return {static_cast<Real>(p[0]), static_cast<Real>(p[1]), static_cast<Real>(p[2])}; }

// An impulse along a row: kept in the row, and given to its bodies as M^-1 J^T lambda. (The velocity solves apply
// the trees' answer this way, from the rows' impulses, not from the bodies' own part of it: the two agree only to
// single precision, which over a chain of light links was a few parts in ten thousand, and what is kept for the
// warm start must be what was applied. This way momentum is kept exactly, too.)
void apply_impulse(JointRow& r, Real lambda) {
    r.impulse += lambda;
    r.a->vel += r.lin_a * (r.a->inv_mass * lambda);
    r.a->w += r.inv_i_ang_a * lambda;
    r.b->vel += r.lin_b * (r.b->inv_mass * lambda);
    r.b->w += r.inv_i_ang_b * lambda;
}

// A rigid row that holds both ways: what the trees solve.
bool two_way(const JointRow& r) { return r.lower == -kInf && r.upper == kInf && r.softness == 0 && r.mass > 0; }
// A one-sided row (a limit, the cone) past its limit, as the position pass builds it: an error below zero.
bool past_limit(const JointRow& r) { return r.lower == 0 && r.upper == kInf && r.error < 0 && r.softness == 0 && r.mass > 0; }

// Rows whose speeds (or errors) answer impulses lambda on them by W lambda, W n x n and symmetric but for rounding:
// lambda to bring them as near as they can come to `need`. First a Cholesky factoring of W, dropping any row whose
// pivot vanishes: that row is implied, to first order, by the rows before it and the tree's (a straight chain fixed
// at both ends cannot be pulled along its length: its cut row along the chain is then such a row), and it is given
// no impulse. That is only to find them. What is solved is N lambda = W need over the rows kept, N = W W (taking W's
// rows and columns of those): lambda by least squares (step 31). When they can all have what they need (as for speeds
// they always can) that is the same lambda as W lambda = need; when they cannot (positions where no pose meets them
// all: a body held where the rows of two joints, both fixed at their other ends, put it at two different places)
// solving W lambda = need over the rows kept met those exactly and left whatever followed for the dropped one, which
// could be much more than all of them had needed: nearly three times as far apart as a pose could bring them.
// factor_least_squares leaves N's factor in L (W made exactly symmetric); solve_least_squares then solves.
void factor_least_squares(int n, std::vector<double>& W, std::vector<double>& L, std::uint8_t* dropped) {
    auto at = [n](int i, int j) { return static_cast<size_t>(i * n + j); };
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j) W[at(i, j)] = W[at(j, i)] = 0.5 * (W[at(i, j)] + W[at(j, i)]);
    L.assign(static_cast<size_t>(n * n), 0.0);
    double scale = 0;
    for (int i = 0; i < n; ++i) scale = std::max(scale, W[at(i, i)]);
    for (int c = 0; c < n; ++c) {
        double d = W[at(c, c)];
        for (int k = 0; k < c; ++k) d -= L[at(c, k)] * L[at(c, k)];
        dropped[c] = !(d > 1e-6 * scale);
        if (dropped[c]) {
            L[at(c, c)] = 1;
            continue;  // (its column below stays zero)
        }
        L[at(c, c)] = std::sqrt(d);
        for (int i = c + 1; i < n; ++i) {
            double e = W[at(i, c)];
            for (int k = 0; k < c; ++k) e -= L[at(i, k)] * L[at(c, k)];
            L[at(i, c)] = e / L[at(c, c)];
        }
    }
    std::vector<double> N(static_cast<size_t>(n * n), 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (dropped[i] || dropped[j]) continue;
            double sum = 0;
            for (int k = 0; k < n; ++k) sum += W[at(i, k)] * W[at(k, j)];
            N[at(i, j)] = sum;
        }
    L.assign(static_cast<size_t>(n * n), 0.0);
    for (int c = 0; c < n; ++c) {
        if (dropped[c]) {
            L[at(c, c)] = 1;
            continue;
        }
        double d = N[at(c, c)];
        for (int k = 0; k < c; ++k) d -= L[at(c, k)] * L[at(c, k)];
        L[at(c, c)] = std::sqrt(std::max(d, 1e-300));
        for (int i = c + 1; i < n; ++i) {
            double e = N[at(i, c)];
            for (int k = 0; k < c; ++k) e -= L[at(i, k)] * L[at(c, k)];
            L[at(i, c)] = e / L[at(c, c)];
        }
    }
}
void solve_least_squares(int n, const std::vector<double>& W, const std::vector<double>& L, const std::uint8_t* dropped,
                         const double* need, double* lambda) {
    auto at = [n](int i, int j) { return static_cast<size_t>(i * n + j); };
    for (int i = 0; i < n; ++i) {  // N lambda = W need, by L L^T (a dropped row gets none)
        lambda[i] = 0;
        if (!dropped[i])
            for (int k = 0; k < n; ++k) lambda[i] += W[at(i, k)] * need[k];
    }
    for (int i = 0; i < n; ++i) {
        if (dropped[i]) {
            lambda[i] = 0;
            continue;
        }
        for (int k = 0; k < i; ++k) lambda[i] -= L[at(i, k)] * lambda[k];
        lambda[i] /= L[at(i, i)];
    }
    for (int i = n - 1; i >= 0; --i) {
        if (dropped[i]) continue;
        for (int k = i + 1; k < n; ++k) lambda[i] -= L[at(k, i)] * lambda[k];
        lambda[i] /= L[at(i, i)];
    }
}

}  // namespace

void JointTrees::child_block(int c, double out[6][6]) const {
    const Node& child = nodes_[static_cast<size_t>(c)];
    const int p = child.parent;
    if (!child.body) {  // a joint below a body: its rows as they act on that body
        const Links& l = links_[static_cast<size_t>(c)];
        const int k = l.node[0] == p ? 0 : 1;
        for (int r = 0; r < child.size; ++r)
            for (int d = 0; d < 6; ++d) out[r][d] = l.J[k][r][d];
    } else {  // a body below a joint: the joint's rows on it, transposed
        const Links& l = links_[static_cast<size_t>(p)];
        const int k = l.node[0] == c ? 0 : 1;
        for (int d = 0; d < 6; ++d)
            for (int r = 0; r < nodes_[static_cast<size_t>(p)].size; ++r) out[d][r] = l.J[k][r][d];
    }
}

int JointTrees::add_node(int size, Body* body) {
    if (used_ == static_cast<int>(nodes_.size())) {
        nodes_.emplace_back();
        links_.emplace_back();
    }
    Node& node = nodes_[static_cast<size_t>(used_)];
    node.size = size;
    node.parent = -1;
    node.body = body;
    links_[static_cast<size_t>(used_)].node[0] = links_[static_cast<size_t>(used_)].node[1] = -1;
    return used_++;
}

void JointTrees::build(const std::vector<Body>& bodies, const std::vector<JointRow>& rows, const std::vector<int>& first,
                       std::vector<std::uint8_t>* taken, bool violated) {
    if (taken) taken->assign(first.empty() ? 0 : first.size() - 1, 0);
    first_ = first;
    base_ = bodies.data();
    // A joint taken in is marked 1, or 2 if one of the rows taken was a limit past its limit.
    auto mark = [&](int joint, int row) {
        std::uint8_t& t = (*taken)[static_cast<size_t>(joint)];
        t = std::max(t, static_cast<std::uint8_t>(past_limit(rows[static_cast<size_t>(row)]) ? 2 : 1));
    };
    used_ = 0;
    steps_.clear();
    data_.clear();
    trees_ = 0;
    nodes_.reserve(first.size() + bodies.size());  // (at most a node per joint and one per body)
    links_.reserve(first.size() + bodies.size());
    const Body* base = bodies.data();
    // A body the trees can take: one that moves, with a mass and inertia that can be inverted.
    auto body_index = [&](const Body* b) -> int {
        const std::ptrdiff_t i = b - base;
        return i >= 0 && static_cast<size_t>(i) < bodies.size() && b->type == BodyType::Dynamic ? static_cast<int>(i) : -1;
    };
    auto usable = [](const Body& b) {
        return b.inv_mass > 0 && b.inv_inertia_body.cx.x > 0 && b.inv_inertia_body.cy.y > 0 && b.inv_inertia_body.cz.z > 0;
    };

    // The nodes: one per joint with two-way rows, one per moving body they hold; except a joint that would close a
    // loop, which is set aside as a cut. Taken in turn, a joint joins two sets of bodies (the world being one more
    // body, so that everything fixed to it is one set); if they are already one, it closes a loop. Joints between
    // two moving bodies are taken first, and joints to the world or a fixed body after: so a structure fixed at both
    // ends is one tree fixed at one, cut at the other (taken in any order, both ends could be fixed first, and the
    // cut would fall in the middle, between two trees, where it cannot be put back through either).
    body_node_.assign(bodies.size(), -1);
    const int world = static_cast<int>(bodies.size());
    union_.resize(bodies.size() + 1);
    for (size_t i = 0; i < union_.size(); ++i) union_[i] = static_cast<int>(i);
    auto find = [&](int i) {
        while (union_[static_cast<size_t>(i)] != i) i = union_[static_cast<size_t>(i)] = union_[static_cast<size_t>(union_[static_cast<size_t>(i)])];
        return i;
    };
    cuts_.clear();
    for (int fixed_side = 0; fixed_side < 2; ++fixed_side)
    for (size_t k = 0; k + 1 < first.size(); ++k) {
        int chosen[6], m = 0;
        for (int r = first[k]; r < first[k + 1]; ++r)
            if ((two_way(rows[static_cast<size_t>(r)]) || (violated && past_limit(rows[static_cast<size_t>(r)]))) && m < 6) chosen[m++] = r;
        if (m == 0) continue;
        const JointRow& r0 = rows[static_cast<size_t>(first[k])];
        const int side[2] = {body_index(r0.a), body_index(r0.b)};
        if ((side[0] < 0 || side[1] < 0) != (fixed_side == 1)) continue;  // (this joint's turn is the other round)
        bool ok = side[0] >= 0 || side[1] >= 0;
        for (int s = 0; s < 2; ++s) ok = ok && (side[s] < 0 || usable(bodies[static_cast<size_t>(side[s])]));
        if (!ok) continue;
        const int set_a = find(side[0] < 0 ? world : side[0]), set_b = find(side[1] < 0 ? world : side[1]);
        if (set_a == set_b) {
            Cut cut;
            cut.joint = static_cast<int>(k);
            cut.size = m;
            std::copy(chosen, chosen + m, cut.rows);
            cut.body[0] = side[0], cut.body[1] = side[1];
            cuts_.push_back(cut);
            continue;
        }
        union_[static_cast<size_t>(set_a)] = set_b;
        const int me = add_node(m, nullptr);
        std::copy(chosen, chosen + m, nodes_[static_cast<size_t>(me)].rows);
        nodes_[static_cast<size_t>(me)].joint = static_cast<int>(k);
        int links = 0;
        for (int s = 0; s < 2; ++s) {
            if (side[s] < 0) continue;
            int& bn = body_node_[static_cast<size_t>(side[s])];
            if (bn < 0) bn = add_node(6, const_cast<Body*>(&bodies[static_cast<size_t>(side[s])]));
            Links& l = links_[static_cast<size_t>(me)];
            l.node[links] = bn;
            for (int i = 0; i < m; ++i) {
                const JointRow& r = rows[static_cast<size_t>(chosen[i])];
                const Vec3 lin = s == 0 ? r.lin_a : r.lin_b, ang = s == 0 ? r.ang_a : r.ang_b;
                const double row[6] = {lin.x, lin.y, lin.z, ang.x, ang.y, ang.z};
                std::copy(row, row + 6, l.J[links][i]);
            }
            ++links;
        }
    }
    if (used_ == 0) return;

    // Who touches whom, as one array: node i's neighbours are adjacency_[adjacent_start_[i] .. [i + 1]).
    const int n = used_;
    adjacent_start_.assign(static_cast<size_t>(n) + 1, 0);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 2; ++k)
            if (const int b = links_[static_cast<size_t>(i)].node[k]; !nodes_[static_cast<size_t>(i)].body && b >= 0)
                ++adjacent_start_[static_cast<size_t>(i) + 1], ++adjacent_start_[static_cast<size_t>(b) + 1];
    for (int i = 0; i < n; ++i) adjacent_start_[static_cast<size_t>(i) + 1] += adjacent_start_[static_cast<size_t>(i)];
    adjacency_.resize(static_cast<size_t>(adjacent_start_[static_cast<size_t>(n)]));
    std::vector<int>& fill = post_;  // (borrowed: the next free slot of each node)
    fill.assign(adjacent_start_.begin(), adjacent_start_.end() - 1);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 2; ++k)
            if (const int b = links_[static_cast<size_t>(i)].node[k]; !nodes_[static_cast<size_t>(i)].body && b >= 0) {
                adjacency_[static_cast<size_t>(fill[static_cast<size_t>(i)]++)] = b;
                adjacency_[static_cast<size_t>(fill[static_cast<size_t>(b)]++)] = i;
            }
    auto neighbours = [&](int i) {
        return std::span<const int>(adjacency_.data() + adjacent_start_[static_cast<size_t>(i)],
                                    adjacency_.data() + adjacent_start_[static_cast<size_t>(i) + 1]);
    };

    // Each connected piece: a tree if it has one edge fewer than nodes. A joint held to the world becomes the
    // root, eliminated last, so that it is never left with only its own (empty) diagonal block. A second one
    // would be (a loop through the world), and the factoring below refuses it.
    seen_.assign(static_cast<size_t>(n), 0);
    std::vector<int>& piece = piece_;
    std::vector<int>& post = post_;
    for (int start = 0; start < n; ++start) {
        if (seen_[static_cast<size_t>(start)]) continue;
        // (the piece, found breadth first: it is its own work list)
        piece.assign(1, start);
        seen_[static_cast<size_t>(start)] = 1;
        for (size_t next = 0; next < piece.size(); ++next)
            for (const int j : neighbours(piece[next]))
                if (!seen_[static_cast<size_t>(j)]) seen_[static_cast<size_t>(j)] = 1, piece.push_back(j);
        size_t edges = 0;
        int held = -1, a_body = -1;
        for (const int i : piece) {
            const Node& node = nodes_[static_cast<size_t>(i)];
            if (node.body) {
                if (a_body < 0) a_body = i;
                continue;
            }
            const int links = (links_[static_cast<size_t>(i)].node[0] >= 0) + (links_[static_cast<size_t>(i)].node[1] >= 0);
            edges += static_cast<size_t>(links);
            if (links == 1) held = i;
        }
        if (edges + 1 != piece.size()) continue;  // (cannot happen: loops were cut above)

        // Leaves first: a depth-first walk from the root, each node listed after everything below it.
        const int root = held >= 0 ? held : a_body;
        post.clear();
        nodes_[static_cast<size_t>(root)].parent = -1;
        walk_.assign(1, {root, 0});
        while (!walk_.empty()) {
            const int i = walk_.back().first;
            const auto around = neighbours(i);
            if (walk_.back().second < around.size()) {
                const int j = around[walk_.back().second++];
                if (j == nodes_[static_cast<size_t>(i)].parent) continue;
                nodes_[static_cast<size_t>(j)].parent = i;
                walk_.push_back({j, 0});
            } else {
                post.push_back(i);
                walk_.pop_back();
            }
        }

        // Factor: each node's diagonal block, less what eliminating its children took from it.
        for (const int i : post) {
            Node& node = nodes_[static_cast<size_t>(i)];
            for (auto& line : node.Dinv) std::fill(std::begin(line), std::end(line), 0.0);
            if (node.body) {
                const double m = 1.0 / node.body->inv_mass;
                const Mat3 I = node.body->inertia_world();
                for (int d = 0; d < 3; ++d) node.Dinv[d][d] = m;
                const Vec3 col[3] = {I.cx, I.cy, I.cz};
                for (int c = 0; c < 3; ++c) {
                    node.Dinv[3][3 + c] = col[c].x;
                    node.Dinv[4][3 + c] = col[c].y;
                    node.Dinv[5][3 + c] = col[c].z;
                }
            }
        }
        bool ok = true;
        for (const int c : post) {
            Node& child = nodes_[static_cast<size_t>(c)];
            if (!kInvert[static_cast<size_t>(child.size)](child.Dinv, child.body ? 1 : -1)) {
                ok = false;
                break;
            }
            if (child.parent < 0) continue;
            Node& parent = nodes_[static_cast<size_t>(child.parent)];
            double B[6][6];
            child_block(c, B);
            kEliminate[static_cast<size_t>(child.size)][static_cast<size_t>(parent.size)](child.Dinv, B, child.T, parent.Dinv);
        }
        if (!ok) continue;
        ++trees_;
        if (taken)
            for (const int i : post)
                if (const Node& node = nodes_[static_cast<size_t>(i)]; !node.body)
                    for (int r = 0; r < node.size; ++r) mark(node.joint, node.rows[r]);
        // Packed for solving, in the same order.
        step_of_.resize(static_cast<size_t>(n));
        for (size_t k = 0; k < post.size(); ++k) step_of_[static_cast<size_t>(post[k])] = static_cast<int>(steps_.size() + k);
        for (const int i : post) {
            const Node& node = nodes_[static_cast<size_t>(i)];
            Step st;
            st.size = node.size;
            st.parent = node.parent < 0 ? -1 : step_of_[static_cast<size_t>(node.parent)];
            st.parent_size = node.parent < 0 ? 0 : nodes_[static_cast<size_t>(node.parent)].size;
            st.data = static_cast<int>(data_.size());
            st.vec = 0;
            st.body = node.body;
            std::copy(node.rows, node.rows + 6, st.rows);
            st.tree = trees_ - 1;  // (counted above)
            st.joint = node.body ? -1 : node.joint;
            for (int a = 0; a < st.size; ++a)
                for (int b = 0; b < st.size; ++b) data_.push_back(node.Dinv[a][b]);
            for (int a = 0; a < st.size; ++a)
                for (int b = 0; b < st.parent_size; ++b) data_.push_back(node.T[a][b]);
            steps_.push_back(st);
        }
    }
    // Interleaved, tree with tree: each tree's nodes are a chain of dependencies (a parent waits for its
    // children), and taken one tree at a time the processor can never work on two nodes at once. Sorted by
    // height instead (leaves 0, a parent one above its highest child), children still come before parents, and
    // the nodes of a level, from every tree, are independent of each other.
    std::vector<int>& height = piece_;  // (borrowed)
    height.assign(steps_.size(), 0);
    for (size_t k = 0; k < steps_.size(); ++k)
        if (const int p = steps_[k].parent; p >= 0) height[static_cast<size_t>(p)] = std::max(height[static_cast<size_t>(p)], height[k] + 1);
    std::vector<int>& by_height = post_;  // (borrowed)
    by_height.resize(steps_.size());
    for (size_t k = 0; k < steps_.size(); ++k) by_height[k] = static_cast<int>(k);
    std::stable_sort(by_height.begin(), by_height.end(), [&](int a, int b) { return height[static_cast<size_t>(a)] < height[static_cast<size_t>(b)]; });
    std::vector<int>& place = step_of_;  // (borrowed: old index -> new)
    place.resize(steps_.size());
    for (size_t k = 0; k < by_height.size(); ++k) place[static_cast<size_t>(by_height[k])] = static_cast<int>(k);
    sorted_.resize(steps_.size());
    for (size_t k = 0; k < by_height.size(); ++k) {
        Step st = steps_[static_cast<size_t>(by_height[k])];
        if (st.parent >= 0) st.parent = place[static_cast<size_t>(st.parent)];
        sorted_[k] = st;
    }
    steps_.swap(sorted_);
    int at = 0;
    for (Step& st : steps_) st.vec = at, at += 2 * st.size;
    vec_.assign(static_cast<size_t>(at), 0.0f);

    // Each tree's own steps, in the same order, and which step each body is.
    tree_start_.assign(static_cast<size_t>(trees_) + 1, 0);
    for (const Step& st : steps_) ++tree_start_[static_cast<size_t>(st.tree) + 1];
    for (int t = 0; t < trees_; ++t) tree_start_[static_cast<size_t>(t) + 1] += tree_start_[static_cast<size_t>(t)];
    tree_steps_.resize(steps_.size());
    std::vector<int>& cursor = piece_;  // (borrowed)
    cursor.assign(tree_start_.begin(), tree_start_.end() - 1);
    body_step_.assign(bodies.size(), -1);
    one_sided_.assign(static_cast<size_t>(trees_), 0);
    for (size_t k = 0; k < steps_.size(); ++k) {
        tree_steps_[static_cast<size_t>(cursor[static_cast<size_t>(steps_[k].tree)]++)] = static_cast<int>(k);
        if (steps_[k].body) body_step_[static_cast<size_t>(steps_[k].body - base)] = static_cast<int>(k);
        if (const int j = steps_[k].joint; j >= 0)
            for (int r = first[static_cast<size_t>(j)]; r < first[static_cast<size_t>(j) + 1]; ++r) {
                const JointRow& row = rows[static_cast<size_t>(r)];
                if (row.lower == 0 && row.upper == kInf && row.softness == 0 && row.mass > 0) one_sided_[static_cast<size_t>(steps_[k].tree)] = 1;
            }
    }

    // The loops: each tree's cut rows, and how its bodies answer an impulse on each.
    loops_.clear();
    response_.clear();
    std::vector<int>& loop_of = post_;  // (borrowed: per tree, its loop, -1 for none, -2 for too many cut rows)
    loop_of.assign(static_cast<size_t>(trees_), -1);
    for (const Cut& cut : cuts_) {
        int step[2] = {-1, -1}, tree = -1;
        bool ok = true;
        for (int s = 0; s < 2; ++s) {
            if (cut.body[s] < 0) continue;
            step[s] = body_step_[static_cast<size_t>(cut.body[s])];
            if (step[s] < 0) ok = false;  // (its tree could not be factored)
            else if (tree < 0) tree = steps_[static_cast<size_t>(step[s])].tree;
            else if (steps_[static_cast<size_t>(step[s])].tree != tree) ok = false;  // (a loop through the world joining two trees)
        }
        if (!ok || tree < 0) continue;
        int& l = loop_of[static_cast<size_t>(tree)];
        if (l == -2) continue;
        if (l == -1) {
            l = static_cast<int>(loops_.size());
            loops_.emplace_back();
            loops_.back().tree = tree;
        }
        Loop& loop = loops_[static_cast<size_t>(l)];
        if (loop.size + cut.size > kMaxCutRows) {
            loop.size = 0;  // (left, all of it, to Gauss-Seidel)
            l = -2;
            continue;
        }
        for (int i = 0; i < cut.size; ++i) {
            const JointRow& r = rows[static_cast<size_t>(cut.rows[i])];
            const int at_row = loop.size++;
            loop.rows[at_row] = cut.rows[i];
            loop.joint[at_row] = cut.joint;
            for (int s = 0; s < 2; ++s) {
                loop.step[at_row][s] = step[s];
                const Vec3 lin = s == 0 ? r.lin_a : r.lin_b, ang = s == 0 ? r.ang_a : r.ang_b;
                const Scalar j[6] = {lin.x, lin.y, lin.z, ang.x, ang.y, ang.z};
                std::copy(j, j + 6, loop.J[at_row][s]);
            }
        }
    }
    for (Loop& loop : loops_) {
        if (loop.size == 0) continue;
        const int first_step = tree_start_[static_cast<size_t>(loop.tree)], end_step = tree_start_[static_cast<size_t>(loop.tree) + 1];
        loop.length = 0;
        for (int k = first_step; k < end_step; ++k) loop.length += static_cast<size_t>(steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])].size);
        loop.response = response_.size();
        response_.resize(response_.size() + loop.length * static_cast<size_t>(loop.size));
        // Column j: the tree solved with cut row j's unit impulse on its bodies (and nothing asked of the tree's
        // own rows): how every node of the tree answers it.
        for (int j = 0; j < loop.size; ++j) {
            for (int k = first_step; k < end_step; ++k) {
                const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
                std::fill(vec_.begin() + st.vec, vec_.begin() + st.vec + st.size, 0.0f);
            }
            for (int s = 0; s < 2; ++s)
                if (loop.step[j][s] >= 0)
                    for (int d = 0; d < 6; ++d) vec_[static_cast<size_t>(steps_[static_cast<size_t>(loop.step[j][s])].vec + d)] += loop.J[j][s][d];
            solve_tree(loop.tree);
            Scalar* column = response_.data() + loop.response + loop.length * static_cast<size_t>(j);
            for (int k = first_step; k < end_step; ++k) {
                const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
                column = std::copy(vec_.begin() + st.vec + st.size, vec_.begin() + st.vec + 2 * st.size, column);
            }
        }
        // W = J_cut A, and its factoring (above). (Its entries: row i of the cut rows, against the bodies' answer to
        // column j, read back from the column at the place of the body's step.)
        std::vector<int>& start_of = step_of_;  // (borrowed: where each of the tree's steps starts in a column)
        start_of.assign(steps_.size(), 0);
        for (int k = first_step, offset = 0; k < end_step; ++k) {
            start_of[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])] = offset;
            offset += steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])].size;
        }
        const int count = loop.size;
        loop.L.assign(static_cast<size_t>(count * count), 0.0);
        std::vector<double>& W = loop.W;
        W.assign(static_cast<size_t>(count * count), 0.0);
        for (int j = 0; j < count; ++j) {
            const Scalar* column = response_.data() + loop.response + loop.length * static_cast<size_t>(j);
            for (int i = 0; i < count; ++i) {
                double sum = 0;
                for (int s = 0; s < 2; ++s)
                    if (loop.step[i][s] >= 0)
                        for (int d = 0; d < 6; ++d) sum += static_cast<double>(loop.J[i][s][d]) * column[start_of[static_cast<size_t>(loop.step[i][s])] + d];
                W[static_cast<size_t>(i * count + j)] = sum;
            }
        }
        factor_least_squares(count, W, loop.L, loop.dropped);
    }
    if (taken)
        for (const Loop& loop : loops_)
            for (int i = 0; i < loop.size; ++i) mark(loop.joint[i], loop.rows[i]);
}

namespace {

// The two halves of a solve for one node, with its size N and its parent's size P fixed at compile time, so that
// the small products are unrolled (there are only a few dozen combinations; see kernels below).
template <int N, int P>
void down(const JointTrees::Scalar* Dinv, const JointTrees::Scalar* y, JointTrees::Scalar* x, JointTrees::Scalar* parent_y) {
    using Scalar = JointTrees::Scalar;
    for (int i = 0; i < N; ++i) {
        Scalar z = 0;
        for (int k = 0; k < N; ++k) z += Dinv[i * N + k] * y[k];
        x[i] = z;
    }
    const Scalar* T = Dinv + N * N;
    for (int j = 0; j < P; ++j) {
        Scalar sum = 0;
        for (int k = 0; k < N; ++k) sum += T[k * P + j] * y[k];
        parent_y[j] -= sum;
    }
}
template <int N, int P>
void up(const JointTrees::Scalar* Dinv, const JointTrees::Scalar* parent_x, JointTrees::Scalar* x) {
    using Scalar = JointTrees::Scalar;
    const Scalar* T = Dinv + N * N;
    for (int i = 0; i < N; ++i) {
        Scalar sum = 0;
        for (int j = 0; j < P; ++j) sum += T[i * P + j] * parent_x[j];
        x[i] -= sum;
    }
}
using Down = void (*)(const JointTrees::Scalar*, const JointTrees::Scalar*, JointTrees::Scalar*, JointTrees::Scalar*);
using Up = void (*)(const JointTrees::Scalar*, const JointTrees::Scalar*, JointTrees::Scalar*);
template <int N, int... P>
constexpr std::array<Down, 7> downs(std::integer_sequence<int, P...>) { return {down<N, P>...}; }
template <int N, int... P>
constexpr std::array<Up, 7> ups(std::integer_sequence<int, P...>) { return {up<N, P>...}; }
using Sizes = std::make_integer_sequence<int, 7>;
// kernels[N][P]: N from 1 to 6, P from 0 (a root) to 6.
constexpr std::array<std::array<Down, 7>, 7> kDown = {downs<1>(Sizes{}), downs<1>(Sizes{}), downs<2>(Sizes{}), downs<3>(Sizes{}),
                                                      downs<4>(Sizes{}), downs<5>(Sizes{}), downs<6>(Sizes{})};
constexpr std::array<std::array<Up, 7>, 7> kUp = {ups<1>(Sizes{}), ups<1>(Sizes{}), ups<2>(Sizes{}), ups<3>(Sizes{}),
                                                  ups<4>(Sizes{}), ups<5>(Sizes{}), ups<6>(Sizes{})};

}  // namespace

void JointTrees::solve_tree(int tree) {
    Scalar* v = vec_.data();
    Scalar none[6];
    const int first_step = tree_start_[static_cast<size_t>(tree)], end_step = tree_start_[static_cast<size_t>(tree) + 1];
    for (int k = first_step; k < end_step; ++k) {
        const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
        kDown[static_cast<size_t>(st.size)][static_cast<size_t>(st.parent_size)](
            data_.data() + st.data, v + st.vec, v + st.vec + st.size, st.parent < 0 ? none : v + steps_[static_cast<size_t>(st.parent)].vec);
    }
    for (int k = end_step - 1; k >= first_step; --k) {
        const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
        if (st.parent < 0) continue;
        const Step& parent = steps_[static_cast<size_t>(st.parent)];
        kUp[static_cast<size_t>(st.size)][static_cast<size_t>(st.parent_size)](data_.data() + st.data, v + parent.vec + parent.size, v + st.vec + st.size);
    }
}

void JointTrees::close_loops(std::vector<JointRow>& rows, bool positions, Real max_linear, Real max_angular) {
    for (const Loop& loop : loops_) {
        const int n = loop.size;
        if (n == 0) continue;
        // What the cut rows still need, after the tree's own correction: for velocities, measured from the bodies
        // (already moved); for positions, their errors less what the tree's correction (still in x) did to them.
        double lambda[kMaxCutRows];
        for (int i = 0; i < n; ++i) {
            const JointRow& row = rows[static_cast<size_t>(loop.rows[i])];
            if (positions) {
                const Real limit = row.lin_a.length_sq() == 0 && row.lin_b.length_sq() == 0 ? max_angular : max_linear;
                double need = -static_cast<double>(std::clamp(row.error, -limit, limit));
                for (int s = 0; s < 2; ++s)
                    if (loop.step[i][s] >= 0) {
                        const Step& st = steps_[static_cast<size_t>(loop.step[i][s])];
                        for (int d = 0; d < 6; ++d) need -= static_cast<double>(loop.J[i][s][d]) * vec_[static_cast<size_t>(st.vec + 6 + d)];
                    }
                lambda[i] = need;
                Real& worst = tree_error_[static_cast<size_t>(loop.tree)];
                worst = std::max(worst, std::fabs(row.error));
            } else {
                lambda[i] = -static_cast<double>(row.velocity() + row.bias);
            }
        }
        double need[kMaxCutRows];
        std::copy(lambda, lambda + n, need);
        solve_least_squares(n, loop.W, loop.L, loop.dropped, need, lambda);
        // A lambda, node by node.
        const int first_step = tree_start_[static_cast<size_t>(loop.tree)], end_step = tree_start_[static_cast<size_t>(loop.tree) + 1];
        size_t offset = 0;
        for (int k = first_step; k < end_step; ++k) {
            const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
            Scalar change[6] = {};
            for (int j = 0; j < n; ++j) {
                const Scalar* column = response_.data() + loop.response + loop.length * static_cast<size_t>(j) + offset;
                for (int d = 0; d < st.size; ++d) change[d] += column[d] * lambda[j];
            }
            offset += static_cast<size_t>(st.size);
            if (st.body && positions) {
                st.body->pos += as_vec3(change);
                st.body->q = integrate_orientation(st.body->q, as_vec3(change + 3), 1);
            } else if (!st.body && !positions) {
                for (int r = 0; r < st.size; ++r) apply_impulse(rows[static_cast<size_t>(st.rows[r])], -static_cast<Real>(change[r]));
            }
        }
        if (!positions)
            for (int i = 0; i < n; ++i) apply_impulse(rows[static_cast<size_t>(loop.rows[i])], static_cast<Real>(lambda[i]));
    }
}

// How far short of its target a contact with a fixed body may be left: no more is solved. In the exact solve (m/s),
// most are left a few millimetres a second short, which the next step's push-out takes care of; a ragdoll's limbs on
// steps were up to 10 m/s short as they landed. In the position pass (m), a tenth of a millimetre.
constexpr double kTouchSlack = 0.01, kTouchSlackPositions = 1e-4;

void JointTrees::correct(std::vector<JointRow>& rows, Real max_linear, Real max_angular, std::vector<Real>& largest,
                         std::span<FixedContact> contacts) {
    if (!contacts.empty()) {
        touch_pos_.resize(body_step_.size());
        touch_q_.resize(body_step_.size());
        for (const FixedContact& c : contacts)
            if (c.body >= 0 && static_cast<size_t>(c.body) < body_step_.size() && body_step_[static_cast<size_t>(c.body)] >= 0) {
                const Body* body = steps_[static_cast<size_t>(body_step_[static_cast<size_t>(c.body)])].body;
                touch_pos_[static_cast<size_t>(c.body)] = body->pos;
                touch_q_[static_cast<size_t>(c.body)] = body->q;
            }
    }
    pass(rows, true, max_linear, max_angular, &largest);
    if (!contacts.empty()) touch(rows, contacts, /*positions=*/true);
}

bool JointTrees::touch(std::vector<JointRow>& rows, std::span<FixedContact> contacts, bool positions) {
    by_tree_.clear();
    for (size_t i = 0; i < contacts.size(); ++i) {
        const int body = contacts[i].body;
        if (body < 0 || static_cast<size_t>(body) >= body_step_.size()) continue;
        const int step = body_step_[static_cast<size_t>(body)];
        if (step >= 0) by_tree_.push_back({steps_[static_cast<size_t>(step)].tree, static_cast<int>(i)});
    }
    std::sort(by_tree_.begin(), by_tree_.end());
    // Every tree, for velocities (its limits and cones may be short with no contact at all); for poses, those with
    // contacts.
    bool changed = false;
    size_t first = 0;
    for (int tree = 0; tree < trees_; ++tree) {
        size_t end = first;
        while (end < by_tree_.size() && by_tree_[end].first == tree) ++end;
        if (end > first || (!positions && one_sided_[static_cast<size_t>(tree)])) changed = touch_tree(rows, contacts, tree, first, end, positions) || changed;
        first = end;
    }
    return changed;
}

bool JointTrees::touch_tree(std::vector<JointRow>& rows, std::span<FixedContact> contacts, int tree, size_t first, size_t end,
                            bool positions) {
    // Where each contact stands now: its speed, or how far the correction has moved it. Nothing to do unless one is
    // short.
    const double slack = positions ? kTouchSlackPositions : kTouchSlack;
    touching_.clear();
    bool any = false;
    for (size_t k = first; k < end; ++k) {
        const FixedContact& c = contacts[static_cast<size_t>(by_tree_[k].second)];
        const int step = body_step_[static_cast<size_t>(c.body)];
        const Body* body = steps_[static_cast<size_t>(step)].body;
        double now;
        if (positions) {
            Quat turn = body->q * conjugate(touch_q_[static_cast<size_t>(c.body)]);
            if (turn.w < 0) turn = turn * Real(-1);
            now = static_cast<double>(dot(c.lin, body->pos - touch_pos_[static_cast<size_t>(c.body)]) + dot(c.ang, turn.vector_part() * 2));
        } else {
            now = static_cast<double>(dot(c.lin, body->vel) + dot(c.ang, body->w));
        }
        const double short_by = static_cast<double>(c.target) - now;
        any = any || short_by > slack;
        touching_.push_back({by_tree_[k].second, -1, {step, -1}, {{c.lin.x, c.lin.y, c.lin.z, c.ang.x, c.ang.y, c.ang.z}, {}}, short_by});
    }
    const int first_step = tree_start_[static_cast<size_t>(tree)], end_step = tree_start_[static_cast<size_t>(tree) + 1];
    // For velocities, the tree's joints' limits and cones too: they may only push as well, and they came before (in
    // the sweeps), so the exact solve could undo them, as it did the contacts (a hinge locked by its limits, joining
    // two bodies flung apart at up to 20 m/s, was left turning at up to 66 rad/s), and so could the contacts. Left out, a ragdoll's forearm and thigh, resting on
    // the steps against their joints' limits, were turned into the limits by the contacts every step and put back by
    // the position pass: turning at 0.04 rad/s without moving, and the ragdoll never fell asleep.
    if (!positions)
        for (int k = first_step; k < end_step; ++k) {
            const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
            if (st.body || st.joint < 0) continue;
            for (int r = first_[static_cast<size_t>(st.joint)]; r < first_[static_cast<size_t>(st.joint) + 1]; ++r) {
                const JointRow& row = rows[static_cast<size_t>(r)];
                if (!(row.lower == 0 && row.upper == kInf && row.softness == 0 && row.mass > 0)) continue;
                Touching t{-1, r, {-1, -1}, {}, -static_cast<double>(row.bias + row.velocity())};
                any = any || t.short_by > slack;
                const Body* side[2] = {row.a, row.b};
                for (int s = 0; s < 2; ++s) {
                    const std::ptrdiff_t i = side[s] - base_;
                    if (i >= 0 && static_cast<size_t>(i) < body_step_.size() && body_step_[static_cast<size_t>(i)] >= 0) {
                        t.step[s] = body_step_[static_cast<size_t>(i)];
                        const Vec3 lin = s == 0 ? row.lin_a : row.lin_b, ang = s == 0 ? row.ang_a : row.ang_b;
                        const double J[6] = {lin.x, lin.y, lin.z, ang.x, ang.y, ang.z};
                        std::copy(J, J + 6, t.J[s]);
                    }
                }
                touching_.push_back(t);
            }
        }
    if (!any) return false;
    // The tree's cut rows come first, if it has a loop: solved again with the contacts, so that they still hold.
    const Loop* loop = nullptr;
    for (const Loop& l : loops_)
        if (l.tree == tree && l.size > 0) loop = &l;
    const int cut = loop ? loop->size : 0;
    const size_t room = static_cast<size_t>(kMaxCutRows - cut);
    if (touching_.size() > room) {  // (the furthest short, if there are too many)
        std::partial_sort(touching_.begin(), touching_.begin() + static_cast<std::ptrdiff_t>(room), touching_.end(),
                          [](const Touching& x, const Touching& y) { return x.short_by > y.short_by; });
        touching_.resize(room);
    }
    const int count = static_cast<int>(touching_.size()), n = cut + count;

    // How the tree's nodes answer a unit impulse on each of them, as build() finds for cut rows.
    touch_offset_.resize(steps_.size());
    size_t length = 0;
    for (int k = first_step; k < end_step; ++k) {
        const int st = tree_steps_[static_cast<size_t>(k)];
        touch_offset_[static_cast<size_t>(st)] = static_cast<int>(length);
        length += static_cast<size_t>(steps_[static_cast<size_t>(st)].size);
    }
    touch_response_.resize(length * static_cast<size_t>(count));
    for (int j = 0; j < count; ++j) {
        for (int k = first_step; k < end_step; ++k) {
            const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
            std::fill(vec_.begin() + st.vec, vec_.begin() + st.vec + st.size, 0.0);
        }
        const Touching& t = touching_[static_cast<size_t>(j)];
        for (int s = 0; s < 2; ++s)
            if (t.step[s] >= 0)
                for (int d = 0; d < 6; ++d) vec_[static_cast<size_t>(steps_[static_cast<size_t>(t.step[s])].vec + d)] += t.J[s][d];
        solve_tree(tree);
        Scalar* column = touch_response_.data() + length * static_cast<size_t>(j);
        for (int k = first_step; k < end_step; ++k) {
            const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
            column = std::copy(vec_.begin() + st.vec + st.size, vec_.begin() + st.vec + 2 * st.size, column);
        }
    }
    auto column_of = [&](int j) -> const Scalar* {
        return j < cut ? response_.data() + loop->response + loop->length * static_cast<size_t>(j)
                       : touch_response_.data() + length * static_cast<size_t>(j - cut);
    };
    // Row i (cut row i, or the other rows' i - cut) against an answer: how much that answer changes its speed.
    auto speed_change = [&](int i, const Scalar* column) {
        double sum = 0;
        for (int s = 0; s < 2; ++s) {
            const int st = i < cut ? loop->step[i][s] : touching_[static_cast<size_t>(i - cut)].step[s];
            if (st < 0) continue;
            const double* J = i < cut ? loop->J[i][s] : touching_[static_cast<size_t>(i - cut)].J[s];
            for (int d = 0; d < 6; ++d) sum += J[d] * column[touch_offset_[static_cast<size_t>(st)] + d];
        }
        return sum;
    };
    std::vector<double> W(static_cast<size_t>(n * n)), need(static_cast<size_t>(n)), lambda(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) W[static_cast<size_t>(i * n + j)] = speed_change(i, column_of(j));
    for (int i = 0; i < cut; ++i) {  // (for poses: left as the correction left them)
        const JointRow& row = rows[static_cast<size_t>(loop->rows[i])];
        need[static_cast<size_t>(i)] = positions ? 0.0 : -static_cast<double>(row.velocity() + row.bias);
    }
    for (int i = cut; i < n; ++i) need[static_cast<size_t>(i)] = touching_[static_cast<size_t>(i - cut)].short_by;
    auto at = [n](int i, int j) { return static_cast<size_t>(i * n + j); };

    // The joints first: the cut rows hold whatever the contacts do. Given the contacts' impulses lambda_S, the cut
    // rows' are lambda_C = W_CC^-1 (need_C - W_CS lambda_S) (by least squares, as close_loops solves them); so the
    // contacts answer their own impulses through A = W_SS - W_SC W_CC^-1 W_CS, and need b = need_S - W_SC W_CC^-1
    // need_C. (Solved as one least-squares system with the contacts, as at first, a contact that could not be met
    // took the joint with it: a body hinged to the world and lying on the floor was left opening its hinge.)
    std::vector<double> through(static_cast<size_t>(cut * (count + 1)));  // W_CC^-1 need_C, then W_CC^-1 W_C(contact j)
    if (cut > 0) {
        double rhs[kMaxCutRows], out[kMaxCutRows];
        for (int k = 0; k <= count; ++k) {
            for (int i = 0; i < cut; ++i) rhs[i] = k == 0 ? need[static_cast<size_t>(i)] : W[at(i, cut + k - 1)];
            solve_least_squares(cut, loop->W, loop->L, loop->dropped, rhs, out);
            std::copy(out, out + cut, through.begin() + static_cast<std::ptrdiff_t>(k * cut));
        }
    }
    std::vector<double> A(static_cast<size_t>(count * count)), b(static_cast<size_t>(count)), low(static_cast<size_t>(count)), x(static_cast<size_t>(count), 0.0);
    std::vector<std::uint8_t> fixed(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i) {
        b[static_cast<size_t>(i)] = need[static_cast<size_t>(cut + i)];
        for (int c = 0; c < cut; ++c) b[static_cast<size_t>(i)] -= W[at(cut + i, c)] * through[static_cast<size_t>(c)];
        for (int j = 0; j < count; ++j) {
            double e = W[at(cut + i, cut + j)];
            for (int c = 0; c < cut; ++c) e -= W[at(cut + i, c)] * through[static_cast<size_t>((1 + j) * cut + c)];
            A[static_cast<size_t>(i * count + j)] = e;
        }
        const Touching& t = touching_[static_cast<size_t>(i)];
        // How the row would answer an impulse with its bodies free, J M^-1 J^T, to compare its answer with the joints
        // against. (Against the most any row of the tree answers, as at first, a contact's was always a thousandth of a
        // forearm's limit's: an impulse turns a light limb far more, in radians a second, than it moves it in metres a
        // second. Every contact of a ragdoll with a limit beside it was left out, by up to 8 m/s short.)
        double free_answer = 0;
        for (int s = 0; s < 2; ++s) {
            if (t.step[s] < 0) continue;
            const Body* body = steps_[static_cast<size_t>(t.step[s])].body;
            const Vec3 lin{static_cast<Real>(t.J[s][0]), static_cast<Real>(t.J[s][1]), static_cast<Real>(t.J[s][2])};
            const Vec3 ang{static_cast<Real>(t.J[s][3]), static_cast<Real>(t.J[s][4]), static_cast<Real>(t.J[s][5])};
            free_answer += static_cast<double>(body->inv_mass * dot(lin, lin) + dot(ang, body->inv_inertia_world() * ang));
        }
        fixed[static_cast<size_t>(i)] = !(A[static_cast<size_t>(i * count + i)] > 1e-5 * free_answer);
        const Real* had = t.contact >= 0 ? contacts[static_cast<size_t>(t.contact)].impulse : &rows[static_cast<size_t>(t.row)].impulse;
        low[static_cast<size_t>(i)] = positions || !had ? 0.0 : -static_cast<double>(*had);  // (in all, never a pull)
    }
    // Then the contacts: each pushes no less than nothing in all, and is either met or not pushing (a contact that
    // is no longer short gives back what it had over). Projected Gauss-Seidel on this small system, from the bodies'
    // exact answer to each: one contact at a time, its impulse what it needs given the others', clamped; until the
    // changes are negligible. (Choosing which are met by trying, adding the furthest short or dropping one that would
    // pull, went round in circles where neither choice was right, the attempts ran out, and a contact was left
    // pulling.) Where the contacts cannot all be met (a body hinged to the world and lying on the floor can only turn
    // one way, and two contacts may want it turned opposite ways), two of them pushing against each other through the
    // hinge could each push harder every sweep without end, impulses in the millions that cancel, and applied in single
    // precision they opened the hinge. So the sweeps stop as soon as they stop helping: when the contacts' shortfall,
    // all of them together, has not fallen by a hundredth in three sweeps. A contact the joints all but fix (its answer
    // to an impulse a hundred thousandth of what it would be with its bodies free, or less) is left as it is: it would
    // take an impulse far out of proportion to move it at all. (A door hinged beside its post, touching it near the
    // hinge, answers a ten thousandth as much as it would free: that one is taken in.)
    double shortfall_was = std::numeric_limits<double>::infinity();
    int stalled = 0;
    for (int sweep = 0; sweep < 200; ++sweep) {
        double biggest = 0;
        for (int i = 0; i < count; ++i) {
            const double d = A[static_cast<size_t>(i * count + i)];
            if (fixed[static_cast<size_t>(i)]) continue;
            double r = b[static_cast<size_t>(i)];
            for (int j = 0; j < count; ++j) r -= A[static_cast<size_t>(i * count + j)] * x[static_cast<size_t>(j)];
            const double next = std::max(low[static_cast<size_t>(i)], x[static_cast<size_t>(i)] + r / d);
            biggest = std::max(biggest, std::fabs(next - x[static_cast<size_t>(i)]) * d);
            x[static_cast<size_t>(i)] = next;
        }
        if (biggest < 1e-7) break;
        double shortfall = 0;
        for (int i = 0; i < count; ++i) {
            if (fixed[static_cast<size_t>(i)]) continue;
            double r = b[static_cast<size_t>(i)];
            for (int j = 0; j < count; ++j) r -= A[static_cast<size_t>(i * count + j)] * x[static_cast<size_t>(j)];
            shortfall += std::max(r, 0.0);
        }
        stalled = shortfall > 0.99 * shortfall_was ? stalled + 1 : 0;
        if (stalled >= 3) break;
        shortfall_was = std::min(shortfall_was, shortfall);
    }
    for (int i = 0; i < count; ++i) lambda[static_cast<size_t>(cut + i)] = x[static_cast<size_t>(i)];
    for (int c = 0; c < cut; ++c) {
        double l = through[static_cast<size_t>(c)];
        for (int j = 0; j < count; ++j) l -= through[static_cast<size_t>((1 + j) * cut + c)] * x[static_cast<size_t>(j)];
        lambda[static_cast<size_t>(c)] = l;
    }

    // Apply it. Poses: each body moved and turned by its part of the answer. Velocities: through the tree's own rows,
    // the cut rows and the contacts, each as an impulse.
    for (int k = first_step; k < end_step; ++k) {
        const Step& st = steps_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
        if (positions != (st.body != nullptr)) continue;
        const int offset = touch_offset_[static_cast<size_t>(tree_steps_[static_cast<size_t>(k)])];
        Scalar change[6] = {};
        for (int r = 0; r < st.size; ++r)
            for (int j = 0; j < n; ++j) change[r] += column_of(j)[offset + r] * lambda[static_cast<size_t>(j)];
        if (positions) {
            st.body->pos += as_vec3(change);
            st.body->q = integrate_orientation(st.body->q, as_vec3(change + 3), 1);
        } else {
            for (int r = 0; r < st.size; ++r) apply_impulse(rows[static_cast<size_t>(st.rows[r])], -static_cast<Real>(change[r]));
        }
    }
    if (positions) return true;
    for (int i = 0; i < cut; ++i) apply_impulse(rows[static_cast<size_t>(loop->rows[i])], static_cast<Real>(lambda[static_cast<size_t>(i)]));
    for (int i = cut; i < n; ++i) {
        const Real l = static_cast<Real>(lambda[static_cast<size_t>(i)]);
        if (l == 0) continue;
        const Touching& t = touching_[static_cast<size_t>(i - cut)];
        if (t.row >= 0) {
            apply_impulse(rows[static_cast<size_t>(t.row)], l);
            continue;
        }
        const FixedContact& c = contacts[static_cast<size_t>(t.contact)];
        Body* body = steps_[static_cast<size_t>(t.step[0])].body;
        body->vel += c.lin * (body->inv_mass * l);
        body->w += body->inv_inertia_world() * c.ang * l;
        *c.impulse += l;
    }
    return true;
}

void JointTrees::pass(std::vector<JointRow>& rows, bool positions, Real max_linear, Real max_angular, std::vector<Real>* largest) {
    Scalar* v = vec_.data();
    tree_error_.assign(static_cast<size_t>(trees_), 0);
    // Right-hand sides: nothing for a body; for a row, the change its speed needs (or, for positions, minus its
    // error, limited).
    for (const Step& st : steps_) {
        Scalar* y = v + st.vec;
        if (st.body)
            for (int r = 0; r < 6; ++r) y[r] = 0;
        else
            for (int r = 0; r < st.size; ++r) {
                const JointRow& row = rows[static_cast<size_t>(st.rows[r])];
                if (positions) {
                    const Real limit = row.lin_a.length_sq() == 0 && row.lin_b.length_sq() == 0 ? max_angular : max_linear;
                    y[r] = -std::clamp(row.error, -limit, limit);
                    Real& worst = tree_error_[static_cast<size_t>(st.tree)];
                    worst = std::max(worst, std::fabs(row.error));
                } else {
                    y[r] = -(row.velocity() + row.bias);
                }
            }
    }
    // Down the tree, leaves first: each node's own share x = Dinv y, and what it leaves its parent. (With B the
    // block coupling child to parent, the parent loses B^T Dinv y = T^T y, Dinv being symmetric.) A root passes
    // nothing on: its parent's size is 0, and it writes to a scratch slot.
    Scalar none[6];
    for (const Step& st : steps_)
        kDown[static_cast<size_t>(st.size)][static_cast<size_t>(st.parent_size)](
            data_.data() + st.data, v + st.vec, v + st.vec + st.size, st.parent < 0 ? none : v + steps_[static_cast<size_t>(st.parent)].vec);
    // Back up, root first: each node less what its parent's answer asks of it, Dinv B x_parent = T x_parent.
    for (auto it = steps_.rbegin(); it != steps_.rend(); ++it) {
        const Step& st = *it;
        if (st.parent < 0) continue;
        const Step& parent = steps_[static_cast<size_t>(st.parent)];
        kUp[static_cast<size_t>(st.size)][static_cast<size_t>(st.parent_size)](data_.data() + st.data, v + parent.vec + parent.size, v + st.vec + st.size);
    }
    // The answer: each body's change of velocity, and each row's impulse (minus mu).
    for (const Step& st : steps_) {
        const Scalar* x = v + st.vec + st.size;
        if (st.body && positions) {
            st.body->pos += as_vec3(x);
            st.body->q = integrate_orientation(st.body->q, as_vec3(x + 3), 1);
        } else if (!st.body && !positions) {
            for (int r = 0; r < st.size; ++r) apply_impulse(rows[static_cast<size_t>(st.rows[r])], -static_cast<Real>(x[r]));
        }
    }
    close_loops(rows, positions, max_linear, max_angular);
    if (largest) {
        for (const Step& st : steps_)
            if (st.joint >= 0) (*largest)[static_cast<size_t>(st.joint)] = tree_error_[static_cast<size_t>(st.tree)];
        for (const Loop& loop : loops_)
            for (int i = 0; i < loop.size; ++i) (*largest)[static_cast<size_t>(loop.joint[i])] = tree_error_[static_cast<size_t>(loop.tree)];
    }
}

}  // namespace phys3d
