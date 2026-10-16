#include <phys3d/broadphase.hpp>

#include <algorithm>

namespace phys3d {

namespace {
// The tree is laid out afresh once a quarter of its leaves have been re-inserted. (Measured on settled piles
// of 1000 and 2000 bodies: waiting for all of them costs two to four times the work per step; rebuilding
// after every single re-insertion gains under a tenth more.)
constexpr std::size_t kRebuildAfter = 4;
// ... or once a search costs a tenth more, per pair found and per leaf, than the first search after the last
// rebuild did (step 20). A rebuild costs about a third of one search (median splits, 0.15 ms for 2000 leaves
// against 0.45 for the search), so one bought at a tenth's wear has paid for itself after three quiet steps.
// (Tried on the benchmark's piles, 250 to 2000 bodies: at a tenth they settle with trees 0.94 to 1.05 times
// as dear to search as fresh ones, at a fifth up to 1.17, at a third up to 1.27.)
constexpr double kWorn = 1.1;
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
    //
    // Counting re-insertions is not enough on its own. A pile that comes to rest with a fifth of its leaves
    // re-inserted since the last rebuild keeps that tree for as long as it rests; and the last few
    // re-insertions before it came to rest, which are the ones it keeps, are made where the pile is densest.
    // Settled piles of 2000 bodies were found taking from 1.06 to 2.2 times the search a fresh tree needs,
    // depending only on their histories. So the cost of each search is watched too.
    const bool rebuild = (churn_ * kRebuildAfter >= n || worn_) && n > 0;
    if (rebuild) {
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
    // The cost of a search, per pair it finds and per leaf it has (the two things a good tree still pays
    // for), against what the same measure was just after the last rebuild. The scene changing (a pile
    // forming, say) moves it too; that only costs a rebuild, which sets a new mark.
    const double cost = static_cast<double>(tests_) / static_cast<double>(out.size() + n);
    if (rebuild) fresh_cost_ = cost;
    worn_ = cost > kWorn * fresh_cost_;
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace phys3d
