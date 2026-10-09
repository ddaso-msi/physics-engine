#include <phys/broadphase.hpp>

#include <algorithm>
#include <numeric>

namespace phys {

namespace {

void add_pair(PairList& out, int i, int j, std::span<const std::uint8_t> immovable) {
    if (immovable[static_cast<size_t>(i)] && immovable[static_cast<size_t>(j)]) return;
    out.push_back({std::min(i, j), std::max(i, j)});
}

}  // namespace

PairList brute_force_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable,
                           std::uint64_t* tests) {
    PairList out;
    std::uint64_t count = 0;
    const int n = static_cast<int>(boxes.size());
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            ++count;
            if (overlap(boxes[static_cast<size_t>(i)], boxes[static_cast<size_t>(j)])) add_pair(out, i, j, immovable);
        }
    if (tests) *tests = count;
    return out;  // already sorted: i ascending, then j ascending
}

PairList SweepAndPrune::find_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable) {
    const size_t n = boxes.size();

    // Bodies are only ever appended or truncated from the end, so the survivors are exactly 0..m-1.
    // Drop the indices that no longer exist, then append the new ones at the back (unsorted).
    const size_t known = std::min(order_.size(), n);
    order_.erase(std::remove_if(order_.begin(), order_.end(), [&](int i) { return static_cast<size_t>(i) >= n; }),
                 order_.end());
    for (size_t i = known; i < n; ++i) order_.push_back(static_cast<int>(i));

    auto left = [&](int i) { return boxes[static_cast<size_t>(i)].lo.x; };

    // Insertion sort: linear when nearly sorted. If it starts shifting far more than a coherent frame
    // would (a teleport, many new bodies), give up and use a full sort instead.
    shifts_ = 0;
    const std::uint64_t shift_limit = 8 * n + 64;
    bool give_up = false;
    for (size_t i = 1; i < n && !give_up; ++i) {
        const int v = order_[i];
        const Real key = left(v);
        size_t j = i;
        while (j > 0 && left(order_[j - 1]) > key) {
            order_[j] = order_[j - 1];
            --j;
            if (++shifts_ > shift_limit) {
                give_up = true;
                break;
            }
        }
        order_[j] = v;  // restores the permutation even when we bail out mid-shift
    }
    if (give_up) std::sort(order_.begin(), order_.end(), [&](int a, int b) { return left(a) < left(b); });

    // Sweep. For body i, every later body whose left edge is <= i's right edge overlaps it in x.
    PairList out;
    tests_ = 0;
    for (size_t i = 0; i < n; ++i) {
        const int bi = order_[i];
        const AABB& box_i = boxes[static_cast<size_t>(bi)];
        for (size_t j = i + 1; j < n; ++j) {
            const int bj = order_[j];
            const AABB& box_j = boxes[static_cast<size_t>(bj)];
            ++tests_;  // counts the comparison that ends the scan too, so the figure reflects real work
            if (box_j.lo.x > box_i.hi.x) break;
            if (box_j.lo.y <= box_i.hi.y && box_i.lo.y <= box_j.hi.y) add_pair(out, bi, bj, immovable);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

PairList TreeBroadphase::find_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable) {
    const size_t n = boxes.size();

    // Keep one proxy per body. Bodies come and go only at the end of the array.
    while (proxy_.size() > n) {
        tree_.destroy_proxy(proxy_.back());
        proxy_.pop_back();
    }
    moves_ = 0;
    for (size_t i = 0; i < proxy_.size(); ++i) {
        if (!contains(tree_.box(proxy_[i]), boxes[i])) {  // escaped its padding: re-insert with fresh padding
            tree_.move_proxy(proxy_[i], fatten(boxes[i], margin_));
            ++moves_;
        }
    }
    for (size_t i = proxy_.size(); i < n; ++i) proxy_.push_back(tree_.create_proxy(fatten(boxes[i], margin_), static_cast<int>(i)));

    // Each movable body asks the tree who overlaps it. A pair of two movable bodies is found from the
    // lower index's query (its real box overlaps the other's padded box whenever the real boxes do);
    // a pair with an immovable body can only come from the movable side.
    PairList out;
    const std::uint64_t before = tree_.visits();
    for (int i = 0; i < static_cast<int>(n); ++i) {
        if (immovable[static_cast<size_t>(i)]) continue;
        tree_.query(boxes[static_cast<size_t>(i)], [&](int j) {
            if (j == i) return;
            if (!immovable[static_cast<size_t>(j)] && j < i) return;  // reported from j's query
            add_pair(out, i, j, immovable);
        });
    }
    tests_ = tree_.visits() - before;
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace phys
