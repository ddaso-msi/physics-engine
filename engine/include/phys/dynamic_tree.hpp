#pragma once
// Stage 7: a dynamic bounding-volume hierarchy (an AABB tree that supports insert, remove and move).
//
// Leaves hold one box each (a "proxy"); every internal node holds the box that encloses its two
// children. Asking "which leaves overlap this box?" then skips whole subtrees whose box misses it,
// so a query costs about log(n) instead of n. Inserting picks the sibling that grows the tree's total
// perimeter least, and AVL-style rotations keep the height near log2(n) even if boxes are inserted in
// sorted order (which would otherwise build a linked list).
#include "aabb.hpp"

#include <cstdint>
#include <vector>

namespace phys {

class DynamicTree {
public:
    static constexpr int kNull = -1;

    // Adds a leaf holding `box` and the caller's `user` value; returns the proxy id.
    int create_proxy(const AABB& box, int user);
    void destroy_proxy(int proxy);
    // Replaces the proxy's box (remove and re-insert).
    void move_proxy(int proxy, const AABB& box);

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
};

}  // namespace phys
