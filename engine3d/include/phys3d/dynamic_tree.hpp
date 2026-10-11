#pragma once
// Stage 10, step 6: a dynamic bounding-volume hierarchy in 3D. The algorithm is the 2D engine's
// (phys/dynamic_tree.hpp), over 3D boxes, with surface area in place of perimeter as the cost of a box.
//
// Leaves hold one box each (a "proxy"); every internal node holds the box that encloses its two
// children. Asking "which leaves overlap this box?" then skips whole subtrees whose box misses it,
// so a query costs about log(n) instead of n. Inserting picks the sibling that grows the tree's total
// surface area least, and AVL-style rotations keep the height near log2(n) even if boxes are inserted in
// sorted order (which would otherwise build a linked list).
#include "aabb.hpp"

#include <cstdint>
#include <vector>

namespace phys3d {

class DynamicTree {
public:
    static constexpr int kNull = -1;

    // Adds a leaf holding `box` and the caller's `user` value; returns the proxy id.
    int create_proxy(const AABB& box, int user);
    void destroy_proxy(int proxy);
    // Replaces the proxy's box (remove and re-insert).
    void move_proxy(int proxy, const AABB& box);

    // Throws away the tree's shape and builds a new one over the same leaves (proxy ids and boxes are
    // kept): split the leaves into two halves at the median along the direction they are most spread out,
    // and do the same to each half. A tree grown one insertion at a time is only as good as each insertion
    // could make it in the tree that was there; this one is laid out with every leaf in view.
    void rebuild();

    const AABB& box(int proxy) const { return nodes_[static_cast<size_t>(proxy)].box; }
    int user(int proxy) const { return nodes_[static_cast<size_t>(proxy)].user; }

    // Calls callback(user) for every leaf whose box overlaps `query_box`.
    template <class Callback>
    void query(const AABB& query_box, Callback&& callback) const {
        if (root_ == kNull) return;
        int stack[kMaxStack];
        int top = 0;
        stack[top++] = root_;
        while (top > 0) {
            const Node& node = nodes_[static_cast<size_t>(stack[--top])];
            ++visits_;
            if (!overlap(node.box, query_box)) continue;
            if (node.is_leaf()) {
                callback(node.user);
            } else {
                stack[top++] = node.child1;
                stack[top++] = node.child2;
            }
        }
    }

    // Calls callback(user_a, user_b) once for every pair of leaves whose boxes overlap.
    //
    // Asking the tree "who overlaps this box?" once per leaf starts every question at the root and walks
    // down through the same upper nodes again and again. Asking the tree about ITSELF does the work once:
    // two subtrees whose boxes do not overlap cannot hold an overlapping pair between them, and are
    // dismissed with one comparison however many leaves they hold. The pairs inside a node are the pairs
    // inside each child plus the pairs between the two children; the pairs between two nodes whose boxes
    // overlap are found by splitting the bigger one and trying again.
    template <class Callback>
    void for_each_overlapping_pair(Callback&& callback) const {
        if (root_ == kNull) return;
        std::vector<PairTask>& stack = pair_stack_;
        stack.clear();
        stack.push_back({root_, kNull});
        while (!stack.empty()) {
            const PairTask task = stack.back();
            stack.pop_back();
            const Node& a = nodes_[static_cast<size_t>(task.a)];
            if (task.b == kNull) {
                if (a.is_leaf()) continue;
                stack.push_back({a.child1, kNull});
                stack.push_back({a.child2, kNull});
                stack.push_back({a.child1, a.child2});
                continue;
            }
            const Node& b = nodes_[static_cast<size_t>(task.b)];
            ++visits_;
            if (!overlap(a.box, b.box)) continue;
            if (a.is_leaf() && b.is_leaf()) {
                callback(a.user, b.user);
            } else if (b.is_leaf() || (!a.is_leaf() && a.height >= b.height)) {
                stack.push_back({a.child1, task.b});
                stack.push_back({a.child2, task.b});
            } else {
                stack.push_back({task.a, b.child1});
                stack.push_back({task.a, b.child2});
            }
        }
    }

    // Calls fn(box, height, is_leaf) for every live node (for debug drawing).
    template <class Fn>
    void for_each_node(Fn&& fn) const {
        for (const Node& n : nodes_)
            if (n.height >= 0) fn(n.box, n.height, n.is_leaf());
    }

    int height() const { return root_ == kNull ? 0 : nodes_[static_cast<size_t>(root_)].height; }
    int proxy_count() const { return proxy_count_; }
    // Total nodes visited by query() so far: a deterministic cost measure for tests and benchmarks.
    std::uint64_t visits() const { return visits_; }

    // Checks every structural invariant (links, heights, enclosing boxes, free list). For tests.
    bool validate() const;

private:
    struct Node {
        AABB box;
        int parent = kNull;  // for a free node: the next free node
        int child1 = kNull, child2 = kNull;
        int height = 0;      // 0 for a leaf, -1 for a free node
        int user = -1;
        bool is_leaf() const { return child1 == kNull; }
    };

    // A balanced tree of 2^32 leaves has height ~45; this is generous without being heap-allocated.
    static constexpr int kMaxStack = 256;

    int allocate_node();
    int build_over(int* leaves, int count);
    void free_node(int node);
    void insert_leaf(int leaf);
    void remove_leaf(int leaf);
    int balance(int node);
    void refit_upwards(int node);
    bool validate_node(int node, int parent, int& leaves) const;

    std::vector<Node> nodes_;
    int root_ = kNull;
    int free_list_ = kNull;
    int proxy_count_ = 0;
    mutable std::uint64_t visits_ = 0;
    // Work still to do in for_each_overlapping_pair. b == kNull means "every pair inside a"; otherwise
    // "every pair with one leaf under a and one under b". (Kept between calls to save allocating it.)
    struct PairTask {
        int a, b;
    };
    mutable std::vector<PairTask> pair_stack_;

};

}  // namespace phys3d
