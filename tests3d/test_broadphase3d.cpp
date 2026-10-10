#include "test.hpp"
#include <phys3d/world.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace phys3d;

namespace {

constexpr double kEps = 1e-4;

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Vec3 vec(Real extent) { return {range(-extent, extent), range(-extent, extent), range(-extent, extent)}; }
    Quat rotation() {
        Quat q;
        do q = Quat{range(-1, 1), range(-1, 1), range(-1, 1), range(-1, 1)};
        while (q.length() < 0.2f);
        return q.normalized();
    }
};

AABB random_box(Lcg& rng, Real extent) {
    const Vec3 c{rng.range(0, extent), rng.range(0, extent), rng.range(0, extent)};
    const Vec3 h{rng.range(0.3f, 0.7f), rng.range(0.3f, 0.7f), rng.range(0.3f, 0.7f)};
    return {c - h, c + h};
}

// n boxes at constant density (about one per 8 m^3), so the volume grows with n.
std::vector<AABB> uniform_scene(int n, Lcg& rng) {
    std::vector<AABB> boxes;
    const Real extent = std::cbrt(8.0f * static_cast<Real>(n));
    for (int i = 0; i < n; ++i) boxes.push_back(random_box(rng, extent));
    return boxes;
}

std::vector<std::uint8_t> flags_for(size_t n) {  // every tenth body is immovable
    std::vector<std::uint8_t> f(n, 0);
    for (size_t i = 0; i < n; i += 10) f[i] = 1;
    return f;
}

bool is_superset(const PairList& big, const PairList& small) {
    return std::includes(big.begin(), big.end(), small.begin(), small.end());
}

}  // namespace

// ---- AABB ----------------------------------------------------------------------------------------

TEST(aabb_of_a_sphere_and_of_turned_boxes) {
    const AABB s = compute_aabb(Body::solid_sphere(0.5f, 1, {1, 2, 3}));
    CHECK_NEAR(s.lo.x, 0.5, kEps); CHECK_NEAR(s.hi.y, 2.5, kEps); CHECK_NEAR(s.lo.z, 2.5, kEps);

    const AABB flat = compute_aabb(Body::solid_box({1, 2, 3}, 1, {0, 0, 0}));
    CHECK_NEAR(flat.hi.x, 1, kEps); CHECK_NEAR(flat.hi.y, 2, kEps); CHECK_NEAR(flat.hi.z, 3, kEps);

    // A quarter turn about z swaps the x and y extents.
    const AABB turned = compute_aabb(Body::solid_box({1, 2, 3}, 1, {0, 0, 0}, Quat::from_axis_angle({0, 0, 1}, kPi / 2)));
    CHECK_NEAR(turned.hi.x, 2, kEps); CHECK_NEAR(turned.hi.y, 1, kEps); CHECK_NEAR(turned.hi.z, 3, kEps);

    // 45 degrees about y: the x and z extents become (hx + hz) / sqrt(2).
    const AABB diag = compute_aabb(Body::solid_box({1, 2, 3}, 1, {5, 0, 0}, Quat::from_axis_angle({0, 1, 0}, kPi / 4)));
    CHECK_NEAR(diag.hi.x - 5, 4 / std::sqrt(2.0), kEps);
    CHECK_NEAR(diag.hi.z, 4 / std::sqrt(2.0), kEps);
    CHECK_NEAR(diag.hi.y, 2, kEps);
}

// For any orientation the box must contain all eight corners, and be tight: some corner touches every side.
TEST(aabb_is_tight_around_a_box_in_any_orientation) {
    Lcg rng;
    for (int trial = 0; trial < 300; ++trial) {
        const Vec3 h{rng.range(0.2f, 2), rng.range(0.2f, 2), rng.range(0.2f, 2)};
        const Body b = Body::solid_box(h, 1, rng.vec(5), rng.rotation());
        const AABB box = compute_aabb(b);
        Vec3 lo{1e9f, 1e9f, 1e9f}, hi{-1e9f, -1e9f, -1e9f};
        for (int i = 0; i < 8; ++i) {
            const Vec3 c = apply(b.transform(), {(i & 1) ? h.x : -h.x, (i & 2) ? h.y : -h.y, (i & 4) ? h.z : -h.z});
            for (int k = 0; k < 3; ++k) {
                lo[k] = std::min(lo[k], c[k]);
                hi[k] = std::max(hi[k], c[k]);
            }
        }
        for (int k = 0; k < 3; ++k) {
            CHECK_NEAR(box.lo[k], lo[k], 1e-3);
            CHECK_NEAR(box.hi[k], hi[k], 1e-3);
        }
    }
}

TEST(aabb_operations_in_3d) {
    const AABB a{{0, 0, 0}, {2, 2, 2}}, touching{{2, 0, 0}, {3, 1, 1}}, inside{{0.5f, 0.5f, 0.5f}, {1, 1, 1}};
    CHECK(overlap(a, touching));  // sharing a face counts
    // Apart along exactly one axis, each axis in turn: all three must be checked.
    CHECK(!overlap(a, AABB{{2.1f, 0, 0}, {3, 1, 1}}));
    CHECK(!overlap(a, AABB{{0, 2.1f, 0}, {1, 3, 1}}));
    CHECK(!overlap(a, AABB{{0, 0, 2.1f}, {1, 1, 3}}));
    CHECK(!overlap(a, AABB{{0, 0, -3}, {1, 1, -0.1f}}));
    CHECK(contains(a, inside) && !contains(inside, a));
    CHECK(!contains(a, AABB{{0.5f, 0.5f, 0.5f}, {1, 1, 2.5f}}));  // pokes out along z only
    const AABB m = merge(a, AABB{{-1, 1, 1}, {1, 5, 3}});
    CHECK_NEAR(m.lo.x, -1, kEps); CHECK_NEAR(m.hi.y, 5, kEps); CHECK_NEAR(m.hi.z, 3, kEps); CHECK_NEAR(m.lo.z, 0, kEps);
    const AABB f = fatten(inside, 0.25f);
    CHECK_NEAR(f.lo.z, 0.25, kEps); CHECK_NEAR(f.hi.x, 1.25, kEps);
    CHECK_NEAR(surface_area(AABB{{0, 0, 0}, {1, 2, 3}}), 2 + 6 + 3, kEps);  // half of 22
}

// ---- pairs ---------------------------------------------------------------------------------------

TEST(pairs_of_two_immovable_bodies_are_never_reported_in_3d) {
    const std::vector<AABB> boxes = {{{0, 0, 0}, {2, 2, 2}}, {{1, 1, 1}, {3, 3, 3}}, {{1.5f, 1.5f, 1.5f}, {2.5f, 2.5f, 2.5f}}};
    const std::vector<std::uint8_t> imm = {1, 1, 0};
    const PairList expected = {{0, 2}, {1, 2}};
    TreeBroadphase tree;
    CHECK(brute_force_pairs(boxes, imm) == expected);
    CHECK(tree.find_pairs(boxes, imm) == expected);
}

TEST(tree_stays_balanced_and_valid_in_3d) {
    DynamicTree line;
    for (int i = 0; i < 1000; ++i) line.create_proxy({{static_cast<Real>(i), 0, 0}, {static_cast<Real>(i) + 0.9f, 1, 1}}, i);
    CHECK(line.validate());
    CHECK(line.height() <= 14);  // log2(1000) is about 10; a naive tree would be 999 deep

    Lcg rng;
    rng.s = 99;
    DynamicTree tree;
    std::vector<int> proxies, users;
    std::vector<AABB> boxes;
    int next_user = 0;
    for (int op = 0; op < 4000; ++op) {
        const Real r = rng.next();
        if (proxies.empty() || r < 0.45f) {
            const AABB b = random_box(rng, 20);
            proxies.push_back(tree.create_proxy(b, next_user));
            boxes.push_back(b);
            users.push_back(next_user++);
        } else if (r < 0.75f) {
            const size_t k = static_cast<size_t>(rng.next() * static_cast<Real>(proxies.size())) % proxies.size();
            boxes[k] = random_box(rng, 20);
            tree.move_proxy(proxies[k], boxes[k]);
        } else {
            const size_t k = static_cast<size_t>(rng.next() * static_cast<Real>(proxies.size())) % proxies.size();
            tree.destroy_proxy(proxies[k]);
            proxies.erase(proxies.begin() + static_cast<long>(k));
            boxes.erase(boxes.begin() + static_cast<long>(k));
            users.erase(users.begin() + static_cast<long>(k));
        }
        if (op % 50 == 0) CHECK(tree.validate());
    }
    CHECK(tree.validate());
    CHECK(tree.proxy_count() == static_cast<int>(proxies.size()));

    // A query returns exactly the leaves a linear scan finds.
    for (int q = 0; q < 200; ++q) {
        const AABB query = fatten(random_box(rng, 20), rng.range(0, 3));
        std::vector<int> got, expected;
        tree.query(query, [&](int user) { got.push_back(user); });
        for (size_t k = 0; k < boxes.size(); ++k)
            if (overlap(boxes[k], query)) expected.push_back(users[k]);
        std::sort(got.begin(), got.end());
        std::sort(expected.begin(), expected.end());
        CHECK(got == expected);
    }
}

TEST(tree_broadphase_touches_the_tree_only_when_a_body_escapes_its_padding_in_3d) {
    Lcg rng;
    auto boxes = uniform_scene(100, rng);
    const auto imm = flags_for(boxes.size());
    TreeBroadphase tree(0.1f);
    tree.find_pairs(boxes, imm);
    for (AABB& b : boxes) b = {b.lo + Vec3{0.01f, 0.01f, 0.01f}, b.hi + Vec3{0.01f, 0.01f, 0.01f}};
    tree.find_pairs(boxes, imm);
    CHECK(tree.last_moves() == 0);
    boxes[3] = {boxes[3].lo + Vec3{0, 0, 5}, boxes[3].hi + Vec3{0, 0, 5}};  // one body leaps along z
    tree.find_pairs(boxes, imm);
    CHECK(tree.last_moves() == 1);
}

// Bodies jitter every frame, a third teleport every tenth frame, the count drops at frame 20 and grows
// again at 35. The tree must always report at least every real overlap, and stay structurally valid.
TEST(tree_broadphase_finds_every_overlap_as_bodies_move_appear_and_vanish_in_3d) {
    Lcg rng;
    rng.s = 5;
    std::vector<AABB> boxes = uniform_scene(200, rng);
    TreeBroadphase tree;
    std::size_t total_pairs = 0;
    for (int frame = 0; frame < 60; ++frame) {
        if (frame > 0)
            for (AABB& b : boxes) {
                const Vec3 d = rng.vec(0.05f);
                b = {b.lo + d, b.hi + d};
            }
        if (frame % 10 == 5)
            for (AABB& b : boxes)
                if (rng.next() < 0.3f) b = random_box(rng, 9);
        if (frame == 20) boxes.resize(80);
        if (frame == 35)
            for (int i = 0; i < 150; ++i) boxes.push_back(random_box(rng, 9));
        const auto imm = flags_for(boxes.size());
        const PairList expected = brute_force_pairs(boxes, imm);
        const PairList got = tree.find_pairs(boxes, imm);
        CHECK(is_superset(got, expected));
        CHECK(std::is_sorted(got.begin(), got.end()));
        CHECK(tree.tree().validate());
        CHECK(tree.tree().proxy_count() == static_cast<int>(boxes.size()));
        total_pairs += expected.size();
    }
    CHECK(total_pairs > 1000);  // the scene really had overlaps to find
}

// Quadrupling the number of evenly spread bodies multiplies brute-force work by 16; the tree must do far better.
TEST(tree_broadphase_work_grows_sub_quadratically_in_3d) {
    auto measure = [](int n, std::uint64_t& brute, std::uint64_t& tree_tests) {
        Lcg rng;
        const auto boxes = uniform_scene(n, rng);
        const std::vector<std::uint8_t> imm(boxes.size(), 0);
        brute_force_pairs(boxes, imm, &brute);
        TreeBroadphase tree;
        tree.find_pairs(boxes, imm);
        tree_tests = tree.last_tests();
    };
    std::uint64_t b1, t1, b2, t2;
    measure(500, b1, t1);
    measure(2000, b2, t2);
    CHECK(static_cast<double>(b2) / static_cast<double>(b1) > 15.0);
    CHECK(static_cast<double>(t2) / static_cast<double>(t1) < 10.0);
    CHECK(t2 < b2 / 10);
}

// ---- inside the World ----------------------------------------------------------------------------

namespace {

void build_scene(World& w) {
    Body floor = Body::fixed_box({20, 0.5f, 20}, {0, -0.5f, 0});
    floor.restitution = 0;
    w.add(floor);
    for (int layer = 0; layer < 3; ++layer) {
        const int m = 3 - layer;
        for (int i = 0; i < m; ++i)
            for (int k = 0; k < m; ++k) {
                Body b = Body::solid_box({0.5f, 0.5f, 0.5f}, 1,
                                         {(static_cast<Real>(i) - static_cast<Real>(m - 1) * 0.5f) * 1.02f, 0.5f + static_cast<Real>(layer) * 1.002f,
                                          (static_cast<Real>(k) - static_cast<Real>(m - 1) * 0.5f) * 1.02f});
                b.restitution = 0;
                w.add(b);
            }
    }
    for (int i = 0; i < 4; ++i) {
        Body ball = Body::solid_sphere(0.4f, 1, {-1 + static_cast<Real>(i) * 0.7f, 6 + static_cast<Real>(i) * 1.3f, 0.3f * static_cast<Real>(i)});
        ball.restitution = 0.4f;
        w.add(ball);
    }
}

}  // namespace

// The broad phase must be invisible to the physics: the same scene stepped with either gives bit-for-bit
// the same bodies, because both feed the narrow phase the same pairs in the same order.
TEST(simulation_is_identical_with_either_broadphase_in_3d) {
    World brute, tree;
    brute.broadphase = BroadphaseKind::BruteForce;
    build_scene(brute);
    build_scene(tree);
    CHECK(tree.broadphase == BroadphaseKind::DynamicTree);  // the default
    for (int step = 0; step < 360; ++step) {
        brute.step(1.0f / 120.0f);
        tree.step(1.0f / 120.0f);
    }
    for (size_t i = 0; i < brute.bodies.size(); ++i) {
        const Body &a = brute.bodies[i], &b = tree.bodies[i];
        CHECK(a.pos.x == b.pos.x && a.pos.y == b.pos.y && a.pos.z == b.pos.z);
        CHECK(a.q.x == b.q.x && a.q.w == b.q.w && a.vel.y == b.vel.y && a.w.z == b.w.z);
    }
    CHECK(brute.stats().contacts == tree.stats().contacts);
    CHECK(brute.stats().contacts > 10);
    CHECK(tree.stats().candidate_pairs >= tree.stats().contacts);
    CHECK(tree.broadphase_tree().validate());
    CHECK(tree.broadphase_tree().proxy_count() == static_cast<int>(tree.bodies.size()));
}

TEST(world_copes_with_truncation_under_the_tree_in_3d) {
    World w;
    build_scene(w);
    for (int step = 0; step < 120; ++step) w.step(1.0f / 120.0f);
    w.truncate(1);
    for (int i = 0; i < 12; ++i) w.add(Body::solid_box({0.4f, 0.4f, 0.4f}, 1, {static_cast<Real>(i % 4) * 1.5f - 2, 0.4f + static_cast<Real>(i / 4), 0}));
    for (int step = 0; step < 120; ++step) w.step(1.0f / 120.0f);
    CHECK(w.broadphase_tree().validate());
    CHECK(w.broadphase_tree().proxy_count() == 13);
    CHECK(!w.contacts().empty());
}
