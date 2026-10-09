#include <phys/dynamic_tree.hpp>

#include <algorithm>

namespace phys {

int DynamicTree::allocate_node() {
    if (free_list_ == kNull) {
        nodes_.emplace_back();
        return static_cast<int>(nodes_.size()) - 1;
    }
    const int id = free_list_;
    Node& n = nodes_[static_cast<size_t>(id)];
    free_list_ = n.parent;  // a free node stores the next free node in `parent`
    n = Node{};
    return id;
}

void DynamicTree::free_node(int id) {
    Node& n = nodes_[static_cast<size_t>(id)];
    n = Node{};
    n.parent = free_list_;
    n.height = -1;
    free_list_ = id;
}

int DynamicTree::create_proxy(const AABB& box, int user) {
    const int leaf = allocate_node();
    nodes_[static_cast<size_t>(leaf)].box = box;
    nodes_[static_cast<size_t>(leaf)].user = user;
    nodes_[static_cast<size_t>(leaf)].height = 0;
    insert_leaf(leaf);
    ++proxy_count_;
    return leaf;
}

void DynamicTree::destroy_proxy(int proxy) {
    remove_leaf(proxy);
    free_node(proxy);
    --proxy_count_;
}

void DynamicTree::move_proxy(int proxy, const AABB& box) {
    remove_leaf(proxy);
    nodes_[static_cast<size_t>(proxy)].box = box;
    insert_leaf(proxy);
}

// Recompute the box and height of `node` and every ancestor, rebalancing on the way up.
void DynamicTree::refit_upwards(int node) {
    while (node != kNull) {
        node = balance(node);
        Node& n = nodes_[static_cast<size_t>(node)];
        const Node& c1 = nodes_[static_cast<size_t>(n.child1)];
        const Node& c2 = nodes_[static_cast<size_t>(n.child2)];
        n.height = 1 + std::max(c1.height, c2.height);
        n.box = merge(c1.box, c2.box);
        node = n.parent;
    }
}

void DynamicTree::insert_leaf(int leaf) {
    nodes_[static_cast<size_t>(leaf)].parent = kNull;
    if (root_ == kNull) {
        root_ = leaf;
        return;
    }

    // Walk down from the root choosing, at each internal node, the child that increases the total
    // perimeter least. Stop when making the new leaf a sibling of the current node is cheaper than
    // descending into either child.
    const AABB leaf_box = nodes_[static_cast<size_t>(leaf)].box;
    int index = root_;
    while (!nodes_[static_cast<size_t>(index)].is_leaf()) {
        const Node& n = nodes_[static_cast<size_t>(index)];
        const AABB& box1 = nodes_[static_cast<size_t>(n.child1)].box;
        const AABB& box2 = nodes_[static_cast<size_t>(n.child2)].box;

        const Real area = perimeter(n.box);
        const Real combined_area = perimeter(merge(n.box, leaf_box));
        const Real cost_here = 2 * combined_area;                    // pair the leaf with this node
        const Real inheritance = 2 * (combined_area - area);         // every ancestor grows by this much if we go deeper

        auto descend_cost = [&](int child, const AABB& child_box) {
            const Real merged = perimeter(merge(leaf_box, child_box));
            return nodes_[static_cast<size_t>(child)].is_leaf() ? merged + inheritance
                                                                : (merged - perimeter(child_box)) + inheritance;
        };
        const Real cost1 = descend_cost(n.child1, box1);
        const Real cost2 = descend_cost(n.child2, box2);

        if (cost_here < cost1 && cost_here < cost2) break;
        index = cost1 < cost2 ? n.child1 : n.child2;
    }

    // Splice a new parent in above the chosen sibling.
    const int sibling = index;
    const int old_parent = nodes_[static_cast<size_t>(sibling)].parent;
    const int new_parent = allocate_node();  // may reallocate nodes_: take no references before this
    Node& np = nodes_[static_cast<size_t>(new_parent)];
    np.parent = old_parent;
    np.user = -1;
    np.box = merge(leaf_box, nodes_[static_cast<size_t>(sibling)].box);
    np.height = nodes_[static_cast<size_t>(sibling)].height + 1;
    np.child1 = sibling;
    np.child2 = leaf;
    nodes_[static_cast<size_t>(sibling)].parent = new_parent;
    nodes_[static_cast<size_t>(leaf)].parent = new_parent;
    if (old_parent != kNull) {
        Node& op = nodes_[static_cast<size_t>(old_parent)];
        (op.child1 == sibling ? op.child1 : op.child2) = new_parent;
    } else {
        root_ = new_parent;
    }

    refit_upwards(nodes_[static_cast<size_t>(leaf)].parent);
}

void DynamicTree::remove_leaf(int leaf) {
    if (leaf == root_) {
        root_ = kNull;
        return;
    }
    const int parent = nodes_[static_cast<size_t>(leaf)].parent;
    const int grand = nodes_[static_cast<size_t>(parent)].parent;
    const Node& p = nodes_[static_cast<size_t>(parent)];
    const int sibling = p.child1 == leaf ? p.child2 : p.child1;

    // The sibling takes the parent's place.
    nodes_[static_cast<size_t>(sibling)].parent = grand;
    if (grand != kNull) {
        Node& g = nodes_[static_cast<size_t>(grand)];
        (g.child1 == parent ? g.child1 : g.child2) = sibling;
        free_node(parent);
        refit_upwards(grand);
    } else {
        root_ = sibling;
        free_node(parent);
    }
    nodes_[static_cast<size_t>(leaf)].parent = kNull;
}

// If one child of `a` is more than one level taller than the other, rotate the taller child up so
// the tree stays shallow. Returns the node now at the top of this subtree.
int DynamicTree::balance(int ia) {
    Node& a = nodes_[static_cast<size_t>(ia)];
    if (a.is_leaf() || a.height < 2) return ia;

    const int ib = a.child1, ic = a.child2;
    Node& b = nodes_[static_cast<size_t>(ib)];
    Node& c = nodes_[static_cast<size_t>(ic)];
    const int imbalance = c.height - b.height;

    auto fix_grandparent = [&](int top) {
        // `top` has taken `a`'s place beneath a's old parent.
        if (nodes_[static_cast<size_t>(top)].parent != kNull) {
            Node& gp = nodes_[static_cast<size_t>(nodes_[static_cast<size_t>(top)].parent)];
            (gp.child1 == ia ? gp.child1 : gp.child2) = top;
        } else {
            root_ = top;
        }
    };

    if (imbalance > 1) {  // right side too tall: lift C above A
        const int iF = c.child1, iG = c.child2;
        Node& f = nodes_[static_cast<size_t>(iF)];
        Node& g = nodes_[static_cast<size_t>(iG)];
        c.child1 = ia;
        c.parent = a.parent;
        a.parent = ic;
        fix_grandparent(ic);
        if (f.height > g.height) {  // the shorter grandchild (G) goes under A
            c.child2 = iF;
            a.child2 = iG;
            g.parent = ia;
            a.box = merge(b.box, g.box);
            c.box = merge(a.box, f.box);
            a.height = 1 + std::max(b.height, g.height);
            c.height = 1 + std::max(a.height, f.height);
        } else {
            c.child2 = iG;
            a.child2 = iF;
            f.parent = ia;
            a.box = merge(b.box, f.box);
            c.box = merge(a.box, g.box);
            a.height = 1 + std::max(b.height, f.height);
            c.height = 1 + std::max(a.height, g.height);
        }
        return ic;
    }

    if (imbalance < -1) {  // left side too tall: lift B above A
        const int iD = b.child1, iE = b.child2;
        Node& d = nodes_[static_cast<size_t>(iD)];
        Node& e = nodes_[static_cast<size_t>(iE)];
        b.child1 = ia;
        b.parent = a.parent;
        a.parent = ib;
        fix_grandparent(ib);
        if (d.height > e.height) {
            b.child2 = iD;
            a.child1 = iE;
            e.parent = ia;
            a.box = merge(c.box, e.box);
            b.box = merge(a.box, d.box);
            a.height = 1 + std::max(c.height, e.height);
            b.height = 1 + std::max(a.height, d.height);
        } else {
            b.child2 = iE;
            a.child1 = iD;
            d.parent = ia;
            a.box = merge(c.box, d.box);
            b.box = merge(a.box, e.box);
            a.height = 1 + std::max(c.height, d.height);
            b.height = 1 + std::max(a.height, e.height);
        }
        return ib;
    }
    return ia;
}

bool DynamicTree::validate_node(int id, int parent, int& leaves) const {
    if (id < 0 || id >= static_cast<int>(nodes_.size())) return false;
    const Node& n = nodes_[static_cast<size_t>(id)];
    if (n.parent != parent || n.height < 0) return false;
    if (n.is_leaf()) {
        ++leaves;
        return n.child2 == kNull && n.height == 0 && n.user >= 0;
    }
    if (n.child2 == kNull) return false;
    const Node& c1 = nodes_[static_cast<size_t>(n.child1)];
    const Node& c2 = nodes_[static_cast<size_t>(n.child2)];
    if (n.height != 1 + std::max(c1.height, c2.height)) return false;
    if (!contains(n.box, c1.box) || !contains(n.box, c2.box)) return false;
    return validate_node(n.child1, id, leaves) && validate_node(n.child2, id, leaves);
}

bool DynamicTree::validate() const {
    int leaves = 0;
    if (root_ != kNull && !validate_node(root_, kNull, leaves)) return false;
    if (leaves != proxy_count_) return false;

    int free_count = 0;
    for (int i = free_list_; i != kNull; i = nodes_[static_cast<size_t>(i)].parent) {
        if (nodes_[static_cast<size_t>(i)].height != -1) return false;
        if (++free_count > static_cast<int>(nodes_.size())) return false;  // cycle
    }
    const int live = leaves == 0 ? 0 : 2 * leaves - 1;  // a full binary tree
    return live + free_count == static_cast<int>(nodes_.size());
}

}  // namespace phys
