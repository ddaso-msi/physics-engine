// Where the 3D engine's time goes. Build the release preset for meaningful numbers:
//   cmake --preset release && cmake --build --preset release && ./build/release/bench3d/phys3d_bench
//
// Each scene drops bodies into a closed room, lets them settle for some seconds with sleeping off (so that
// every body is simulated every step), then averages the time per step and per phase over two more seconds.
//   crates   boxes only
//   mixed    spheres, boxes, capsules and hulls in equal numbers
//   lumps    hulls only (twelve random points each)
//   chains   hanging chains of ball-jointed links, swinging into each other
#include <phys3d/world.hpp>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <vector>

using namespace phys3d;

namespace {

struct Lcg {
    std::uint32_t s = 1;
    Real next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<Real>(s >> 8) / static_cast<Real>(1u << 24);
    }
    Real range(Real lo, Real hi) { return lo + (hi - lo) * next(); }
    Quat rotation() { return Quat::from_axis_angle({range(-1, 1), range(-1, 1), range(-1, 1)}, range(-3, 3)); }
};

enum class Kind { Crates, Mixed, Lumps, Chains };

// A floor and four walls `half` from the middle, tall enough that nothing bounces out.
void add_room(World& w, Real half) {
    const Real tall = 200;
    w.add(Body::fixed_box({half + 1, 0.5f, half + 1}, {0, -0.5f, 0}));
    w.add(Body::fixed_box({0.5f, tall, half + 1}, {-half - 0.5f, tall, 0}));
    w.add(Body::fixed_box({0.5f, tall, half + 1}, {half + 0.5f, tall, 0}));
    w.add(Body::fixed_box({half + 1, tall, 0.5f}, {0, tall, -half - 0.5f}));
    w.add(Body::fixed_box({half + 1, tall, 0.5f}, {0, tall, half + 0.5f}));
}

Body make(int which, Lcg& rng, Vec3 at) {
    const Quat q = rng.rotation();
    switch (which) {
        case 0: return Body::solid_sphere(rng.range(0.3f, 0.5f), 1, at);
        case 1: return Body::solid_box({rng.range(0.25f, 0.5f), rng.range(0.25f, 0.5f), rng.range(0.25f, 0.5f)}, 1, at, q);
        case 2: return Body::solid_capsule(rng.range(0.2f, 0.4f), rng.range(0.2f, 0.3f), 1, at, q);
        default: {
            std::vector<Vec3> points;
            for (int k = 0; k < 12; ++k) points.push_back({rng.range(-0.5f, 0.5f), rng.range(-0.5f, 0.5f), rng.range(-0.5f, 0.5f)});
            const std::optional<Body> lump = Body::solid_hull(points.data(), 12, 1, at, q);
            return lump ? *lump : Body::solid_sphere(0.4f, 1, at);
        }
    }
}

World build(Kind kind, int count, NarrowPhase path) {
    World w;
    w.narrowphase = path;
    w.allow_sleep = false;
    Lcg rng;
    if (kind == Kind::Chains) {
        // Chains of 10 links hung in a grid, each started leaning so that neighbours collide.
        const int chains = count / 10, side = static_cast<int>(std::ceil(std::sqrt(static_cast<Real>(chains))));
        add_room(w, static_cast<Real>(side) * 1.2f + 2);
        for (int c = 0; c < chains; ++c) {
            const Vec3 top{(static_cast<Real>(c % side) - static_cast<Real>(side) * 0.5f) * 2.4f, 9, (static_cast<Real>(c / side) - static_cast<Real>(side) * 0.5f) * 2.4f};
            const Vec3 lean = Vec3{rng.range(-1, 1), -1, rng.range(-1, 1)}.normalized();
            const Quat q = Quat::from_axis_angle(cross({0, -1, 0}, lean), std::acos(-lean.y));
            int previous = -1;
            for (int i = 0; i < 10; ++i) {
                const int link = w.add(Body::solid_box({0.08f, 0.35f, 0.08f}, 1, top + lean * (0.35f + 0.7f * static_cast<Real>(i)), q));
                w.add_joint(Joint::ball(w.bodies, previous, link, top + lean * (0.7f * static_cast<Real>(i))));
                previous = link;
            }
        }
        return w;
    }
    // A square room sized so the settled pile is three or four bodies deep, filled from a grid above it.
    const int side = std::max(4, static_cast<int>(std::ceil(std::sqrt(static_cast<Real>(count) / 3.5f))));
    add_room(w, static_cast<Real>(side) * 0.55f + 0.5f);
    for (int i = 0; i < count; ++i) {
        const int layer = i / (side * side), cell = i % (side * side);
        const Vec3 at{(static_cast<Real>(cell % side) - static_cast<Real>(side - 1) * 0.5f) * 1.1f, 1 + static_cast<Real>(layer) * 1.3f,
                      (static_cast<Real>(cell / side) - static_cast<Real>(side - 1) * 0.5f) * 1.1f};
        Body b = make(kind == Kind::Crates ? 1 : kind == Kind::Lumps ? 3 : i % 4, rng, at);
        b.restitution = 0;
        w.add(b);
    }
    return w;
}

void measure(const char* name, Kind kind, int count, NarrowPhase path) {
    World w = build(kind, count, path);
    const Real dt = 1.0f / 120.0f;
    for (int i = 0; i < 6 * 120; ++i) w.step(dt);
    const int frames = 240;
    double total[7] = {};
    std::size_t contacts = 0, points = 0;
    for (int i = 0; i < frames; ++i) {
        w.step(dt);
        const StepStats& s = w.stats();
        const double phase[7] = {s.ms_broadphase, s.ms_narrowphase, s.ms_islands, s.ms_solver, s.ms_movement, s.ms_sweep, s.ms_sleep};
        for (int k = 0; k < 7; ++k) total[k] += phase[k] / frames;
        contacts += s.contacts;
        for (const ContactPair& c : w.contacts()) points += static_cast<std::size_t>(c.manifold.count);
    }
    double step = 0;
    for (const double t : total) step += t;
    std::printf("%-7s %5d %-9s %7.3f | %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f | %6zu %6zu\n", name, count,
                path == NarrowPhase::Convex ? "GJK/EPA" : "per shape", step, total[0], total[1], total[2], total[3], total[4], total[5], total[6],
                contacts / static_cast<std::size_t>(frames), points / static_cast<std::size_t>(frames));
}

}  // namespace

int main(int argc, char** argv) {
    const int largest = argc > 1 ? std::atoi(argv[1]) : 2000;
    std::printf("ms per step at 120 Hz, all bodies awake, 10 solver sweeps\n");
    std::printf("%-7s %5s %-9s %7s | %6s %6s %6s %6s %6s %6s %6s | %6s %6s\n", "scene", "n", "narrow", "step", "broad", "narrow", "island", "solver", "move",
                "sweep", "sleep", "pairs", "points");
    for (int n : {250, 500, 1000, 2000, 4000}) {
        if (n > largest) break;
        measure("crates", Kind::Crates, n, NarrowPhase::Specialised);
        measure("crates", Kind::Crates, n, NarrowPhase::Convex);
        measure("mixed", Kind::Mixed, n, NarrowPhase::Specialised);
        measure("lumps", Kind::Lumps, n, NarrowPhase::Specialised);
        measure("chains", Kind::Chains, n, NarrowPhase::Specialised);
    }
    return 0;
}
