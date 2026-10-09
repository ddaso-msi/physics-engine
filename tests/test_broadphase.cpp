#include "test.hpp"
#include <phys/world.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace phys;

namespace {

constexpr double kEps = 1e-5;

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
};

AABB random_box(Lcg& rng, Real extent) {
    const Real x = rng.range(0, extent), y = rng.range(0, extent);
    const Real w = rng.range(0.3f, 0.7f), h = rng.range(0.3f, 0.7f);
    return {{x - w, y - h}, {x + w, y + h}};
}

// n boxes at constant density (about one per 4 m^2), so the area grows with n.
std::vector<AABB> uniform_scene(int n, Lcg& rng) {
    std::vector<AABB> boxes;
    const Real extent = std::sqrt(4.0f * static_cast<Real>(n));
    for (int i = 0; i < n; ++i) boxes.push_back(random_box(rng, extent));
    return boxes;
}

// Every tenth body is immovable, like scenery.
std::vector<std::uint8_t> flags_for(size_t n) {
    std::vector<std::uint8_t> f(n, 0);
    for (size_t i = 0; i < n; i += 10) f[i] = 1;
    return f;
}

bool is_superset(const PairList& big, const PairList& small) {
    return std::includes(big.begin(), big.end(), small.begin(), small.end());
}

}  // namespace

// ---- AABB ----------------------------------------------------------------------------------------

TEST(aabb_of_a_circle_and_of_rotated_polygons) {
    Body c(Shape::make_circle(0.5f), {1, 2}, 0);
    AABB cb = compute_aabb(c);
    CHECK_NEAR(cb.lo.x, 0.5, kEps); CHECK_NEAR(cb.lo.y, 1.5, kEps);
    CHECK_NEAR(cb.hi.x, 1.5, kEps); CHECK_NEAR(cb.hi.y, 2.5, kEps);

    const Polygon box = Polygon::box(1.0f, 0.5f);
    AABB flat = compute_aabb(Body(Shape::make_polygon(box), {0, 0}, 0));
    CHECK_NEAR(flat.hi.x, 1, kEps); CHECK_NEAR(flat.hi.y, 0.5, kEps);
    AABB upright = compute_aabb(Body(Shape::make_polygon(box), {0, 0}, kPi / 2));
    CHECK_NEAR(upright.hi.x, 0.5, 1e-4); CHECK_NEAR(upright.hi.y, 1, 1e-4);
    AABB diamond = compute_aabb(Body(Shape::make_polygon(Polygon::box(1, 1)), {3, 0}, kPi / 4));
    CHECK_NEAR(diamond.hi.x, 3 + std::sqrt(2.0), 1e-4);  // a unit-half-extent square turned 45 degrees
    CHECK_NEAR(diamond.lo.y, -std::sqrt(2.0), 1e-4);
}

TEST(aabb_operations) {
    const AABB a{{0, 0}, {2, 2}}, touching{{2, 0}, {3, 1}}, apart{{2.1f, 0}, {3, 1}}, inside{{0.5f, 0.5f}, {1, 1}};
    CHECK(overlap(a, touching));  // sharing an edge counts
    CHECK(!overlap(a, apart));
    CHECK(contains(a, inside) && !contains(inside, a));
    const AABB m = merge(a, apart);
    CHECK_NEAR(m.lo.x, 0, kEps); CHECK_NEAR(m.hi.x, 3, kEps); CHECK_NEAR(m.hi.y, 2, kEps);
    const AABB f = fatten(inside, 0.25f);
    CHECK_NEAR(f.lo.x, 0.25, kEps); CHECK_NEAR(f.hi.y, 1.25, kEps);
    CHECK_NEAR(perimeter(a), 4, kEps);
}

// ---- sweep and prune ------------------------------------------------------------------------------

TEST(sweep_and_prune_matches_brute_force) {
    for (int n : {0, 1, 2, 50, 300}) {
        Lcg rng;
        rng.s = static_cast<std::uint32_t>(n) + 7;
        const auto boxes = uniform_scene(n, rng);
        const auto imm = flags_for(boxes.size());
        SweepAndPrune sap;
        CHECK(sap.find_pairs(boxes, imm) == brute_force_pairs(boxes, imm));
    }
}

TEST(pairs_of_two_immovable_bodies_are_never_reported) {
    const std::vector<AABB> boxes = {{{0, 0}, {2, 2}}, {{1, 1}, {3, 3}}, {{1.5f, 1.5f}, {2.5f, 2.5f}}};
    const std::vector<std::uint8_t> imm = {1, 1, 0};  // 0 and 1 overlap but both are scenery
    const PairList expected = {{0, 2}, {1, 2}};
    SweepAndPrune sap;
    TreeBroadphase tree;
    CHECK(brute_force_pairs(boxes, imm) == expected);
    CHECK(sap.find_pairs(boxes, imm) == expected);
    CHECK(tree.find_pairs(boxes, imm) == expected);
}

// ---- dynamic tree ---------------------------------------------------------------------------------

TEST(tree_stays_balanced_when_inserting_in_sorted_order) {
    // A naive tree built from boxes inserted left to right degenerates into a 1000-long chain.
    DynamicTree tree;
    for (int i = 0; i < 1000; ++i)
        tree.create_proxy({{static_cast<Real>(i), 0}, {static_cast<Real>(i) + 0.9f, 1}}, i);
    CHECK(tree.validate());
    CHECK(tree.proxy_count() == 1000);
    CHECK(tree.height() <= 14);  // log2(1000) is about 10
}

TEST(tree_survives_random_inserts_moves_and_removals) {
    Lcg rng;
    rng.s = 99;
    DynamicTree tree;
    std::vector<int> proxies;       // live proxy ids
    std::vector<AABB> boxes;        // their boxes (parallel)
    std::vector<int> users;         // their user values (parallel)
    int next_user = 0;

    for (int op = 0; op < 4000; ++op) {
        const Real r = rng.next();
        if (proxies.empty() || r < 0.45f) {
            const AABB b = random_box(rng, 40);
            proxies.push_back(tree.create_proxy(b, next_user));
            boxes.push_back(b);
            users.push_back(next_user++);
        } else if (r < 0.75f) {
            const size_t k = static_cast<size_t>(rng.next() * static_cast<Real>(proxies.size())) % proxies.size();
            const AABB b = random_box(rng, 40);
            tree.move_proxy(proxies[k], b);
            boxes[k] = b;
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
    CHECK(tree.height() <= 2 * static_cast<int>(std::ceil(std::log2(static_cast<double>(proxies.size()) + 1))) + 2);

    // A query must return exactly the leaves a linear scan finds.
    for (int q = 0; q < 200; ++q) {
        const AABB query = fatten(random_box(rng, 40), rng.range(0, 3));
        std::vector<int> got, expected;
        tree.query(query, [&](int user) { got.push_back(user); });
        for (size_t k = 0; k < boxes.size(); ++k)
            if (overlap(boxes[k], query)) expected.push_back(users[k]);
        std::sort(got.begin(), got.end());
        std::sort(expected.begin(), expected.end());
        CHECK(got == expected);
    }
}

TEST(tree_broadphase_only_touches_the_tree_when_a_body_escapes_its_padding) {
    Lcg rng;
    auto boxes = uniform_scene(100, rng);
    const auto imm = flags_for(boxes.size());
    TreeBroadphase tree(0.1f);
    tree.find_pairs(boxes, imm);

    for (AABB& b : boxes) b = {{b.lo.x + 0.01f, b.lo.y + 0.01f}, {b.hi.x + 0.01f, b.hi.y + 0.01f}};
    tree.find_pairs(boxes, imm);
    CHECK(tree.last_moves() == 0);  // 1 cm of drift fits inside 10 cm of padding

    boxes[3] = {{boxes[3].lo.x + 5, boxes[3].lo.y}, {boxes[3].hi.x + 5, boxes[3].hi.y}};
    tree.find_pairs(boxes, imm);
    CHECK(tree.last_moves() == 1);
}

// ---- all three broad phases, on a changing scene --------------------------------------------------

// Bodies jitter every frame, 30% teleport every tenth frame, the body count drops at frame 20 and
// grows again at 35. Sweep and prune must match brute force exactly; the tree must report at least
// every real overlap.
TEST(broadphases_agree_with_brute_force_as_bodies_move_appear_and_vanish) {
    Lcg rng;
    rng.s = 5;
    std::vector<AABB> boxes = uniform_scene(200, rng);
    SweepAndPrune sap;
    TreeBroadphase tree;

    for (int frame = 0; frame < 60; ++frame) {
        if (frame > 0)
            for (AABB& b : boxes) {
                const Vec2 d{rng.range(-0.05f, 0.05f), rng.range(-0.05f, 0.05f)};
                b = {b.lo + d, b.hi + d};
            }
        if (frame % 10 == 5)
            for (AABB& b : boxes)
                if (rng.next() < 0.3f) b = random_box(rng, 30);
        if (frame == 20) boxes.resize(80);
        if (frame == 35)
            for (int i = 0; i < 150; ++i) boxes.push_back(random_box(rng, 30));

        const auto imm = flags_for(boxes.size());
        const PairList expected = brute_force_pairs(boxes, imm);
        CHECK(sap.find_pairs(boxes, imm) == expected);
        CHECK(is_superset(tree.find_pairs(boxes, imm), expected));
        CHECK(tree.tree().validate());
        CHECK(tree.tree().proxy_count() == static_cast<int>(boxes.size()));
    }
}

// ---- cost ----------------------------------------------------------------------------------------

// Quadrupling the number of evenly spread bodies multiplies brute-force work by 16. A broad phase that
// is worth having must grow much more slowly.
TEST(broadphase_work_grows_sub_quadratically) {
    auto measure = [](int n, std::uint64_t& brute, std::uint64_t& sap_tests, std::uint64_t& tree_tests) {
        Lcg rng;
        const auto boxes = uniform_scene(n, rng);
        const std::vector<std::uint8_t> imm(boxes.size(), 0);
        brute_force_pairs(boxes, imm, &brute);
        SweepAndPrune sap;
        sap.find_pairs(boxes, imm);
        sap_tests = sap.last_tests();
        TreeBroadphase tree;
        tree.find_pairs(boxes, imm);
        tree_tests = tree.last_tests();
    };
    std::uint64_t b1, s1, t1, b2, s2, t2;
    measure(500, b1, s1, t1);
    measure(2000, b2, s2, t2);

    CHECK(static_cast<double>(b2) / static_cast<double>(b1) > 15.0);  // exactly quadratic
    CHECK(static_cast<double>(s2) / static_cast<double>(s1) < 10.0);
    CHECK(static_cast<double>(t2) / static_cast<double>(t1) < 10.0);
    CHECK(s2 < b2 / 10);
    CHECK(t2 < b2 / 10);
}

// Sweep and prune sorts along x only, so a vertical column of boxes (all sharing one x range) makes it
// compare every pair. The tree does not care which way things are arranged.
TEST(tree_beats_sweep_and_prune_on_a_tall_column) {
    const int n = 1000;
    std::vector<AABB> boxes;
    for (int i = 0; i < n; ++i) boxes.push_back({{0, static_cast<Real>(i)}, {1, static_cast<Real>(i) + 0.9f}});
    const std::vector<std::uint8_t> imm(boxes.size(), 0);

    SweepAndPrune sap;
    const PairList sap_pairs = sap.find_pairs(boxes, imm);
    TreeBroadphase tree;
    const PairList tree_pairs = tree.find_pairs(boxes, imm);

    CHECK(sap_pairs.empty() && tree_pairs.size() >= sap_pairs.size());
    CHECK(sap.last_tests() == static_cast<std::uint64_t>(n) * (n - 1) / 2);  // every pair, as predicted
    CHECK(tree.last_tests() < 60000);                                          // about 25 per body
}

// ---- inside the World ------------------------------------------------------------------------------

namespace {

Body crate_at(Vec2 p, Real half = 0.5f) {
    Body b(Shape::make_polygon(Polygon::box(half, half)), p, 0);
    b.restitution = 0;
    return b;
}

void build_test_scene(World& w) {
    Body ground(Shape::make_polygon(Polygon::box(20, 0.5f)), {0, -0.5f}, 0, BodyType::Static);
    ground.restitution = 0;
    w.add(ground);
    for (int row = 0; row < 5; ++row)
        for (int i = 0; i < 5 - row; ++i)
            w.add(crate_at({(static_cast<Real>(i) - static_cast<Real>(4 - row) * 0.5f) * 1.02f,
                            0.5f + static_cast<Real>(row) * 1.002f}));
    for (int i = 0; i < 4; ++i) {
        Body ball(Shape::make_circle(0.4f), {-1.5f + static_cast<Real>(i), 7.0f + static_cast<Real>(i) * 1.3f}, 0);
        ball.restitution = 0.4f;
        w.add(ball);
    }
}

}  // namespace

// The broad phase must be invisible to the physics: the same scene stepped with each of the three
// gives bit-for-bit the same bodies, because they all feed the narrow phase the same pairs in order.
TEST(simulation_is_identical_whichever_broadphase_runs) {
    World brute, sap, tree;
    brute.broadphase = BroadphaseKind::BruteForce;
    sap.broadphase = BroadphaseKind::SweepAndPrune;
    tree.broadphase = BroadphaseKind::DynamicTree;
    build_test_scene(brute);
    build_test_scene(sap);
    build_test_scene(tree);

    for (int step = 0; step < 360; ++step) {
        brute.step(1.0f / 120.0f);
        sap.step(1.0f / 120.0f);
        tree.step(1.0f / 120.0f);
    }
    for (size_t i = 0; i < brute.bodies.size(); ++i) {
        const Body &a = brute.bodies[i], &b = sap.bodies[i], &c = tree.bodies[i];
        CHECK(a.pos.x == b.pos.x && a.pos.y == b.pos.y && a.angle == b.angle && a.vel.x == b.vel.x && a.w == b.w);
        CHECK(a.pos.x == c.pos.x && a.pos.y == c.pos.y && a.angle == c.angle && a.vel.x == c.vel.x && a.w == c.w);
    }
    CHECK(brute.stats().contacts == tree.stats().contacts);
    CHECK(brute.stats().contacts > 10);  // the scene really was in contact, not trivially identical
    // The broad phases earn their keep: far fewer comparisons than brute force.
    CHECK(sap.stats().broadphase_tests < brute.stats().broadphase_tests);
    CHECK(tree.stats().candidate_pairs >= tree.stats().contacts);
}

TEST(world_copes_with_bodies_being_truncated_and_re_added) {
    World w;
    build_test_scene(w);
    for (int step = 0; step < 120; ++step) w.step(1.0f / 120.0f);
    w.truncate(1);
    CHECK(w.bodies.size() == 1);
    for (int i = 0; i < 12; ++i) w.add(crate_at({static_cast<Real>(i) - 6.0f, 0.5f + static_cast<Real>(i % 3)}));
    for (int step = 0; step < 120; ++step) w.step(1.0f / 120.0f);
    CHECK(w.broadphase_tree().validate());
    CHECK(w.broadphase_tree().proxy_count() == 13);
    CHECK(!w.contacts().empty());
}

// ---- collision filtering --------------------------------------------------------------------------

TEST(filter_requires_both_bodies_to_admit_each_other) {
    Body a = crate_at({0, 0}), b = crate_at({0, 0});
    CHECK(should_collide(a, b));  // defaults: everyone collides with everyone

    a.category = 0b10; a.mask = 0b10;
    b.category = 0b01; b.mask = 0b01;
    CHECK(!should_collide(a, b));  // disjoint groups

    b.category = 0b10; b.mask = 0b10;
    CHECK(should_collide(a, b));   // same group

    // One-sided: a refuses b's group while b would accept a's. Still no collision, so a body can
    // never be pushed by something it cannot push back.
    a.mask = 0b01; b.category = 0b10; b.mask = 0xFFFF; a.category = 0b01;
    CHECK(!should_collide(a, b));
}

TEST(filtered_pairs_pass_through_each_other_but_still_hit_the_ground) {
    World w;
    w.gravity = {};
    Body a = crate_at({0, 0}), b = crate_at({0.3f, 0});  // overlapping
    a.category = 0b10; a.mask = 0b10;
    b.category = 0b01; b.mask = 0b01;
    w.add(a);
    w.add(b);
    w.step(1.0f / 60.0f);
    CHECK(w.contacts().empty());
    CHECK(w.stats().candidate_pairs == 1);  // the broad phase still proposed it; the filter dropped it

    World g;
    Body ground(Shape::make_polygon(Polygon::box(10, 0.5f)), {0, -0.5f}, 0, BodyType::Static);
    g.add(ground);
    Body ghost = crate_at({0, 2.0f});
    ghost.category = 0b10;
    ghost.mask = 0b01;  // collides with group 1 (the ground), not with other ghosts
    g.add(ghost);
    for (int step = 0; step < 360; ++step) g.step(1.0f / 120.0f);
    CHECK_NEAR(g.bodies[1].pos.y, 0.5, 0.02);  // landed on the ground
}
