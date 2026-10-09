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
