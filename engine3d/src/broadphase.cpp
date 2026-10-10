#include <phys3d/broadphase.hpp>

#include <algorithm>

namespace phys3d {

namespace {
void add_pair(PairList& out, int i, int j, std::span<const std::uint8_t> immovable) {
    if (immovable[static_cast<size_t>(i)] && immovable[static_cast<size_t>(j)]) return;
    out.push_back({std::min(i, j), std::max(i, j)});
}
}  // namespace

PairList brute_force_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable, std::uint64_t* tests) {
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

PairList TreeBroadphase::find_pairs(std::span<const AABB> boxes, std::span<const std::uint8_t> immovable) {
    const size_t n = boxes.size();

    // One proxy per body; bodies come and go only at the end of the array.
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

    // Each movable body asks the tree who overlaps it. A pair of two movable bodies is reported from the
    // lower index's query; a pair with an immovable body can only come from the movable side.
    PairList out;
    const std::uint64_t before = tree_.visits();
    for (int i = 0; i < static_cast<int>(n); ++i) {
        if (immovable[static_cast<size_t>(i)]) continue;
        tree_.query(boxes[static_cast<size_t>(i)], [&](int j) {
            if (j == i) return;
            if (!immovable[static_cast<size_t>(j)] && j < i) return;
            add_pair(out, i, j, immovable);
        });
    }
    tests_ = tree_.visits() - before;
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace phys3d
