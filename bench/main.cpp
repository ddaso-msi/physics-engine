// Broad phase benchmark. Build the release preset for meaningful numbers:
//   cmake --preset release && cmake --build --preset release && ./build/release/bench/phys_bench
//
// Part 1: broad phase alone, evenly scattered boxes that jitter a little each frame.
// Part 2: the same on a tall column of boxes (sweep and prune's worst case).
// Part 3: whole World::step with crates falling and settling, per broad phase.
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

void part3() {
    std::printf("Whole World::step (ms per step), crates scattered above a floor, 20 steps measured after 60 settling\n");
    std::printf("%7s %10s %12s %10s %10s\n", "bodies", "contacts", "brute", "sweep+prune", "tree");
    for (int n : {100, 250, 500, 1000, 2000}) {
        double ms[3] = {0, 0, 0};
        std::size_t contacts = 0;
        const BroadphaseKind kinds[3] = {BroadphaseKind::BruteForce, BroadphaseKind::SweepAndPrune, BroadphaseKind::DynamicTree};
        for (int k = 0; k < 3; ++k) {
            Lcg rng;
            World w;
            w.broadphase = kinds[k];
            const Real half = std::sqrt(static_cast<Real>(n)) * 0.8f + 2;  // floor half-width
            Body floor(Shape::make_polygon(Polygon::box(half, 0.5f)), {half, -0.5f}, 0, BodyType::Static);
            w.add(floor);
            for (int i = 0; i < n; ++i) {
                Body b(Shape::make_polygon(Polygon::box(0.3f, 0.3f)), {rng.range(0.5f, 2 * half - 0.5f), rng.range(0.5f, 6.0f)}, rng.range(0, 3));
                b.restitution = 0;
                w.add(b);
            }
            for (int i = 0; i < 60; ++i) w.step(1.0f / 60.0f);
            ms[k] = time_ms(0, 20, [&] { w.step(1.0f / 60.0f); });
            contacts = w.stats().contacts;
        }
        std::printf("%7d %10zu %12.3f %10.3f %10.3f\n", n, contacts, ms[0], ms[1], ms[2]);
    }
}

}  // namespace

int main() {
    part1_and_2();
    part3();
    return 0;
}
