// Broad phase benchmark. Build the release preset for meaningful numbers:
//   cmake --preset release && cmake --build --preset release && ./build/release/bench/phys_bench
//
// Part 1: broad phase alone, evenly scattered boxes that jitter a little each frame.
// Part 2: the same on a tall column of boxes (sweep and prune's worst case).
// Part 3: whole World::step with a settled pile of crates, per broad phase.
// Part 4: what sleeping saves on the same pile.
#include <phys/world.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <vector>

using namespace phys;
using Clock = std::chrono::steady_clock;

namespace {

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
};

// Average milliseconds per call over `frames` calls, after `warmup` untimed ones.
double time_ms(int warmup, int frames, const std::function<void()>& fn) {
    for (int i = 0; i < warmup; ++i) fn();
    const auto start = Clock::now();
    for (int i = 0; i < frames; ++i) fn();
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count() / frames;
}

void jitter(std::vector<AABB>& boxes, Lcg& rng) {
    for (AABB& b : boxes) {
        const Vec2 d{rng.range(-0.02f, 0.02f), rng.range(-0.02f, 0.02f)};
        b = {b.lo + d, b.hi + d};
    }
}

void part1_and_2() {
    std::printf("Broad phase only (ms per frame)\n");
    std::printf("%-18s %7s %9s %10s %10s %10s %12s %12s\n", "scene", "bodies", "pairs", "brute", "sweep+prune", "tree",
                "sap tests", "tree visits");
    for (int column = 0; column < 2; ++column) {
        for (int n : {100, 250, 500, 1000, 2000, 5000, 10000}) {
            Lcg rng;
            std::vector<AABB> boxes;
            if (column) {
                for (int i = 0; i < n; ++i) boxes.push_back({{0, static_cast<Real>(i)}, {1, static_cast<Real>(i) + 0.9f}});
            } else {
                const Real extent = std::sqrt(4.0f * static_cast<Real>(n));
                for (int i = 0; i < n; ++i) {
                    const Real x = rng.range(0, extent), y = rng.range(0, extent), w = rng.range(0.3f, 0.7f), h = rng.range(0.3f, 0.7f);
                    boxes.push_back({{x - w, y - h}, {x + w, y + h}});
                }
            }
            const std::vector<std::uint8_t> imm(boxes.size(), 0);
            SweepAndPrune sap;
            TreeBroadphase tree;
            PairList pairs;

            // Each timed call also nudges every box, so SAP's re-sort and the tree's re-insertions are included.
            const double brute_ms = n <= 5000 ? time_ms(1, 5, [&] { jitter(boxes, rng); pairs = brute_force_pairs(boxes, imm); }) : -1;
            const double sap_ms = time_ms(3, 20, [&] { jitter(boxes, rng); pairs = sap.find_pairs(boxes, imm); });
            const double tree_ms = time_ms(3, 20, [&] { jitter(boxes, rng); pairs = tree.find_pairs(boxes, imm); });

            char brute_text[32];
            if (brute_ms < 0) std::snprintf(brute_text, sizeof brute_text, "skipped");
            else std::snprintf(brute_text, sizeof brute_text, "%.3f", brute_ms);
            std::printf("%-18s %7d %9zu %10s %10.3f %10.3f %12llu %12llu\n", column ? "tall column" : "scattered", n,
                        pairs.size(), brute_text, sap_ms, tree_ms, static_cast<unsigned long long>(sap.last_tests()),
                        static_cast<unsigned long long>(tree.last_tests()));
        }
        std::printf("\n");
    }
}

// n unit-ish crates on a jittered grid (so none start overlapping) dropped into a CLOSED room: floor and
// two tall walls. Walls matter: without them crates knocked sideways fall off the world forever, and a
// benchmark that is mostly measuring free fall tells you nothing about contacts.
void build_room(World& w, int n) {
    Lcg rng;
    const Real width = 1.6f * std::sqrt(static_cast<Real>(n)) + 4;
    auto wall = [&](Vec2 half, Vec2 pos) {
        Body b(Shape::make_polygon(Polygon::box(half.x, half.y)), pos, 0, BodyType::Static);
        b.restitution = 0;
        w.add(b);
    };
    wall({width * 0.5f + 1, 0.5f}, {width * 0.5f, -0.5f});
    wall({0.5f, 60}, {-0.5f, 30});
    wall({0.5f, 60}, {width + 0.5f, 30});
    const int cols = static_cast<int>((width - 1) / 0.8f);
    for (int i = 0; i < n; ++i) {
        Body b(Shape::make_polygon(Polygon::box(0.3f, 0.3f)),
               {0.6f + static_cast<Real>(i % cols) * 0.8f, 0.4f + static_cast<Real>(i / cols) * 0.8f}, rng.range(-0.2f, 0.2f));
        b.restitution = 0;
        w.add(b);
    }
}

void part3() {
    std::printf("Whole World::step (ms per step): crates dropped into a closed room, settled for 20 s, sleeping OFF\n");
    std::printf("%7s %10s %12s %12s %10s\n", "bodies", "contacts", "brute", "sweep+prune", "tree");
    for (int n : {250, 500, 1000, 2000}) {
        double ms[3] = {0, 0, 0};
        std::size_t contacts = 0;
        const BroadphaseKind kinds[3] = {BroadphaseKind::BruteForce, BroadphaseKind::SweepAndPrune, BroadphaseKind::DynamicTree};
        for (int k = 0; k < 3; ++k) {
            World w;
            w.allow_sleep = false;
            w.broadphase = kinds[k];
            build_room(w, n);
            for (int i = 0; i < 20 * 60; ++i) w.step(1.0f / 60.0f);
            ms[k] = time_ms(0, 30, [&] { w.step(1.0f / 60.0f); });
            contacts = w.stats().contacts;
        }
        std::printf("%7d %10zu %12.3f %12.3f %10.3f\n", n, contacts, ms[0], ms[1], ms[2]);
    }
    std::printf("\n");
}

void part4() {
    std::printf("Sleeping (tree broad phase): same room, settled for 20 s, then ms per step\n");
    std::printf("%7s %34s %34s\n", "bodies", "sleeping off", "sleeping on");
    for (int n : {250, 500, 1000, 2000}) {
        char text[2][64];
        for (int mode = 0; mode < 2; ++mode) {
            World w;
            w.allow_sleep = mode == 1;
            build_room(w, n);
            for (int i = 0; i < 20 * 60; ++i) w.step(1.0f / 60.0f);
            const double ms = time_ms(0, 30, [&] { w.step(1.0f / 60.0f); });
            std::snprintf(text[mode], sizeof text[mode], "%.3f (%zu awake, %zu islands)", ms, w.stats().awake_bodies, w.stats().islands);
        }
        std::printf("%7d %34s %34s\n", n, text[0], text[1]);
    }
}

}  // namespace

int main() {
    part1_and_2();
    part3();
    part4();
    return 0;
}
