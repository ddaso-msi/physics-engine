#include <phys3d/broadphase.hpp>

#include <algorithm>

namespace phys3d {

namespace {
// The tree is laid out afresh once a quarter of its leaves have been re-inserted. (Measured on settled piles
// of 1000 and 2000 bodies: waiting for all of them costs two to four times the work per step; rebuilding
// after every single re-insertion gains under a tenth more.)
constexpr std::size_t kRebuildAfter = 4;
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
    for (size_t i = proxy_.size(); i < n; ++i) {
        proxy_.push_back(tree_.create_proxy(fatten(boxes[i], margin_), static_cast<int>(i)));
        ++churn_;
    }
    churn_ += static_cast<size_t>(moves_);

    // A tree that has had its leaves taken out and put back one at a time is a worse tree than one laid out
    // fresh over the same boxes: each re-insertion is placed as well as it can be in the tree as it stands,
    // and nothing ever revisits the choices made around it. So once enough leaves have been re-inserted,
    // lay it out again. A pile at rest moves no leaves and so never pays for this.
    if (churn_ * kRebuildAfter >= n && n > 0) {
        tree_.rebuild();
        churn_ = 0;
        ++rebuilds_;
    }

    // Every pair of padded boxes that overlap, from one pass of the tree over itself. The padding admits
    // pairs whose real boxes do not touch, so check those too: what is left is exactly the brute-force list.
    PairList out;
    const std::uint64_t before = tree_.visits();
    tree_.for_each_overlapping_pair([&](int i, int j) {
        if (overlap(boxes[static_cast<size_t>(i)], boxes[static_cast<size_t>(j)])) add_pair(out, i, j, immovable);
    });
    tests_ = tree_.visits() - before;
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace phys3d
