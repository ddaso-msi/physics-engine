#pragma once
// Level 2, Bridge Builder. A car has to cross a gap between two cliffs. The player joins a grid of
// connection points with a limited number of beams; every beam end is a real hinge (its axis across the
// road), to the cliff where the point is fixed in rock, otherwise to the other beams that meet there. Beams
// along the road are the deck the car drives on; the rest are bracing, a single bar under the middle.
//
// A hinge only holds what is pinned to it, and only triangles keep their shape. The right cliff offers no
// fixed point at road level, so a plain row of beams hangs loose at that end: the deck has to be braced back
// to the three points that do hold. Counting freedoms, the smallest bridge that stands is ten beams.
#include "level.hpp"

#include <vector>

namespace puzzle3d {

class Bridge : public Level {
public:
    static constexpr int kColumns = 5, kRows = 2;  // connection points; row 1 is the road
    static constexpr Real kLeftEdge = -4.0f, kSpacing = 2.0f, kRoadY = 3.9f, kLowY = 1.9f;
    static constexpr Real kCliffTop = 4.0f;
    static constexpr Real kDeckHalfWidth = 1.3f;
    static constexpr int kMaxBeams = 12;
    static constexpr Real kFinishX = 8.0f;  // the car's body has to get this far onto the right cliff...
    static constexpr Real kFallY = 2.6f;    // ...without ever dropping this low

    struct Beam {
        int a, b;  // connection point indices, a < b
        bool operator==(const Beam&) const = default;
    };
    struct Car {
        int chassis = -1;
        int wheels[4] = {-1, -1, -1, -1};
        int drives[4] = {-1, -1, -1, -1};  // their hinges, indices into world.joints
    };
    std::vector<Beam> beams;
    Car car;

    Bridge() { reset(); }

    const char* name() const override { return "Bridge Builder"; }
    const char* objective() const override { return "Get the car across the gap."; }
    std::span<const char* const> help() const override;
    const char* hint() const override;
    std::string stats() const override;
    std::string result() const override;
    View view() const override { return {{0, 3.0f, 0}, 0.45f, 0.28f, 17.0f}; }
    Role role(size_t body) const override;
    void overlay(Overlay& out, const Ray* pointer) const override;
    void press(const Ray& ray) override;
    void drag(const Ray& ray) override { pointer_ = ray; }
    void release(const Ray& ray) override;
    bool cancel() override;

    static int node(int column, int row) { return column * kRows + row; }
    static Vec3 node_position(int n);
    static bool is_anchor(int n);
    static bool can_connect(int a, int b);
    static bool is_road(const Beam& b) { return b.a % kRows == 1 && b.b % kRows == 1; }

    bool add_beam(int a, int b);
    bool remove_beam(int a, int b);
    const Body& chassis() const { return world.bodies[static_cast<size_t>(car.chassis)]; }
    // Which connection point the ray passes closest to, if within reach.
    static int node_near(const Ray& ray);

protected:
    void clear_setup() override { beams.clear(); }
    void build() override;
    void evaluate() override;

private:
    void set_throttle(Real speed, Real torque);
    int beam_under(const Ray& ray) const;

    int from_ = -1;
    Ray pointer_{};
    size_t first_beam_ = 0;
};

}  // namespace puzzle3d
