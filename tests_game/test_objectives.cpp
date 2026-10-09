#include "test.hpp"
#include <objectives.hpp>

using namespace phys;
using namespace puzzle;

namespace {
int add_disc(World& w, Vec2 p) { return w.add(Body(Shape::make_circle(0.2f), p, 0)); }
}  // namespace

TEST(count_below_counts_only_listed_bodies_under_the_line) {
    World w;
    const int low = add_disc(w, {0, 1});
    const int high = add_disc(w, {0, 3});
    const int unlisted = add_disc(w, {0, 0});
    (void)unlisted;
    const int listed[] = {low, high};
    CHECK(count_below(w, listed, 2.0f) == 1);
    CHECK(count_below(w, listed, 0.5f) == 0);
    CHECK(count_below(w, listed, 5.0f) == 2);
}

TEST(count_below_uses_the_centre_and_excludes_the_line_itself) {
    World w;
    const int on_line[] = {add_disc(w, {0, 2.0f})};
    CHECK(count_below(w, on_line, 2.0f) == 0);  // radius 0.2 reaches below, the centre does not
    w.bodies[0].pos.y = 1.999f;
    CHECK(count_below(w, on_line, 2.0f) == 1);
}

TEST(count_below_follows_the_simulation) {
    World w;
    const int ball[] = {add_disc(w, {0, 3})};
    CHECK(count_below(w, ball, 2.0f) == 0);
    for (int i = 0; i < 120; ++i) w.step(kTimeStep);  // one second of free fall: 4.9 m
    CHECK(count_below(w, ball, 2.0f) == 1);
}

TEST(inside_tests_the_centre_against_the_region) {
    World w;
    const int a = add_disc(w, {2, 2});
    CHECK(inside(w, a, {1, 1}, {3, 3}));
    CHECK(!inside(w, a, {2.1f, 1}, {3, 3}));
    CHECK(!inside(w, a, {1, 1}, {3, 1.9f}));
    CHECK(inside(w, a, {2, 2}, {2, 2}));  // edges count
}

TEST(all_slow_needs_every_body_under_both_limits) {
    World w;
    const int bodies[] = {add_disc(w, {0, 0}), add_disc(w, {2, 0})};
    CHECK(all_slow(w, bodies));
    w.bodies[1].vel = {0.3f, 0};
    CHECK(!all_slow(w, bodies));
    CHECK(all_slow(w, bodies, 0.5f));
    w.bodies[1].vel = {};
    w.bodies[0].w = 1.0f;  // spinning on the spot is not at rest
    CHECK(!all_slow(w, bodies));
    CHECK(all_slow(w, std::span<const int>{}));
}

TEST(top_of_is_the_highest_point_of_the_outline) {
    World w;
    const int disc = add_disc(w, {0, 3});
    CHECK_NEAR(top_of(w, std::span<const int>{&disc, 1}), 3.2, 1e-5);
    // A 2 x 0.4 plank stood on end reaches a metre above its centre.
    const int plank = w.add(Body(Shape::make_polygon(Polygon::box(1.0f, 0.2f)), {5, 4}, kPi / 2));
    const int both[] = {disc, plank};
    CHECK_NEAR(top_of(w, both), 5.0, 1e-4);
    CHECK_NEAR(top_of(w, std::span<const int>{}), kFloorTop, 1e-6);
}

TEST(hold_needs_an_unbroken_run) {
    Hold h;
    for (int i = 0; i < 59; ++i) CHECK(!h.update(true, kTimeStep, 0.5));
    CHECK(!h.update(false, kTimeStep, 0.5));  // one bad step and the count starts again
    CHECK(h.held == 0);
    int steps = 0;
    while (!h.update(true, kTimeStep, 0.5)) ++steps;
    CHECK(steps >= 58 && steps <= 60);
    CHECK(!h.update(false, kTimeStep, 0.5));
}
