#include "levels/bridge.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace puzzle {

namespace {

constexpr int kCarGroup = 3;
constexpr std::uint16_t kCarBit = 1u << kCarGroup, kDeckBit = 1u << 4, kBraceBit = 1u << 5;
constexpr Real kReach = 0.45f;  // how close the pointer has to be to a connection point

}  // namespace

std::span<const char* const> Bridge::help() const {
    static constexpr const char* kLines[] = {
        "Drag from one connection point to a neighbour to add",
        "a beam. Click a beam to take it away. Twelve at most.",
        "Solid points are fixed in rock. The right cliff's edge",
        "is crumbling: only the point low on its face holds.",
        "Every beam end is a hinge, so only triangles keep their",
        "shape. Beams along the road carry the car.",
    };
    return kLines;
}

const char* Bridge::hint() const {
    if (status == Status::Won) return "Solved. Space to try it with fewer beams.";
    if (status == Status::Failed) return chassis().pos.y < kFallY ? "The car fell into the gap." : "The car did not make it across.";
    if (phase == Phase::Setup) return from_ >= 0 ? "Let go on a neighbouring point." : "Drag between points to add beams, then Space to drive.";
    return "Space to stop and rebuild.";
}

std::string Bridge::stats() const {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "Beams %d/%d", static_cast<int>(beams.size()), kMaxBeams);
    return buf;
}

std::string Bridge::result() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%d beam%s, %.1f s", static_cast<int>(beams.size()), beams.size() == 1 ? "" : "s", win_time);
    return buf;
}

Vec2 Bridge::node_position(int n) {
    return {kLeftEdge + kSpacing * static_cast<Real>(n / kRows), kRoadY + kSpacing * static_cast<Real>(n % kRows - 1)};
}

bool Bridge::is_anchor(int n) {
    const int column = n / kRows, row = n % kRows;
    // The right cliff's edge is crumbling: only the point low on its face holds.
    return (column == 0 && row <= 1) || (column == kColumns - 1 && row == 0);
}

bool Bridge::can_connect(int a, int b) {
    if (a == b || a < 0 || b < 0 || a >= kColumns * kRows || b >= kColumns * kRows) return false;
    if (is_anchor(a) && is_anchor(b)) return false;  // both ends in rock: it would do nothing
    return std::abs(a / kRows - b / kRows) <= 1 && std::abs(a % kRows - b % kRows) <= 1;
}

bool Bridge::add_beam(int a, int b) {
    if (phase != Phase::Setup || !can_connect(a, b) || static_cast<int>(beams.size()) >= kMaxBeams) return false;
    const Beam beam{std::min(a, b), std::max(a, b)};
    if (std::find(beams.begin(), beams.end(), beam) != beams.end()) return false;
    beams.push_back(beam);
    rebuild();
    return true;
}

bool Bridge::remove_beam(int a, int b) {
    if (phase != Phase::Setup) return false;
    const auto it = std::find(beams.begin(), beams.end(), Beam{std::min(a, b), std::max(a, b)});
    if (it == beams.end()) return false;
    beams.erase(it);
    rebuild();
    return true;
}

int Bridge::node_near(Vec2 p) const {
    for (int n = 0; n < kColumns * kRows; ++n)
        if (distance(p, node_position(n)) <= kReach) return n;
    return -1;
}

int Bridge::beam_near(Vec2 p) const {
    for (size_t i = 0; i < beams.size(); ++i)
        if (world.bodies[first_beam_ + i].contains(p)) return static_cast<int>(i);
    return -1;
}

void Bridge::press(Vec2 p) {
    if (phase != Phase::Setup) return;
    from_ = node_near(p);
    if (from_ >= 0) return;
    const int beam = beam_near(p);
    if (beam >= 0) remove_beam(beams[static_cast<size_t>(beam)].a, beams[static_cast<size_t>(beam)].b);
}

void Bridge::release(Vec2 p) {
    if (from_ < 0) return;
    const int from = from_;
    from_ = -1;
    add_beam(from, node_near(p));
}

bool Bridge::cancel() {
    if (from_ < 0) return false;
    from_ = -1;
    return true;
}

Role Bridge::role(size_t body) const {
    if (body >= first_beam_) return Role::Tool;
    if (world.bodies[body].type == phys::BodyType::Static) return Role::Scenery;
    return status == Status::Won ? Role::Done : Role::Load;
}

void Bridge::overlay(Overlay& out, Vec2 pointer, bool pointer_active) const {
    out.lines.push_back({{kFinishX, kCliffTop}, {kFinishX, kCliffTop + 1.6f}, Role::Target, false});
    out.labels.push_back({{kFinishX + 0.15f, kCliffTop + 1.6f}, "finish", Role::Target});
    if (phase != Phase::Setup) return;
    for (int n = 0; n < kColumns * kRows; ++n) {
        const bool lit = n == from_ || (from_ >= 0 && can_connect(from_, n));
        out.dots.push_back({node_position(n), is_anchor(n) ? 0.16f : 0.11f, lit ? Role::Selected : is_anchor(n) ? Role::Scenery : Role::Guide});
    }
    if (from_ >= 0 && pointer_active) out.lines.push_back({node_position(from_), pointer, Role::Selected, true});
}

void Bridge::build() {
    using phys::scenes::add_static_box;
    add_room(world);
    const Real right_edge = kLeftEdge + kSpacing * static_cast<Real>(kColumns - 1);
    const Real h = (kCliffTop - kFloorTop) * 0.5f;
    add_static_box(world, {(kLeftEdge - 0.5f) * 0.5f, h}, {(kLeftEdge + 0.5f) * 0.5f, kFloorTop + h});
    add_static_box(world, {(kWorldW - 0.5f - right_edge) * 0.5f, h}, {(kWorldW - 0.5f + right_edge) * 0.5f, kFloorTop + h});

    car = phys::scenes::build_car(world, {2.6f, kCliffTop + 0.62f + phys::scenes::kCarWheelRadius + 0.01f}, kCarGroup);
    phys::scenes::set_car_throttle(world, car, 7.0f, 25.0f);

    // Beams. Those along the road collide with the car and nothing else; bracing collides with nothing.
    first_beam_ = world.bodies.size();
    for (const Beam& beam : beams) {
        const Vec2 a = node_position(beam.a), b = node_position(beam.b), d = b - a;
        Body body(phys::Shape::make_polygon(phys::Polygon::box(d.length() * 0.5f, 0.1f)), (a + b) * 0.5f, std::atan2(d.y, d.x),
                  phys::BodyType::Dynamic, 2.0f);
        const bool road = beam.a % kRows == 1 && beam.b % kRows == 1;
        body.category = road ? kDeckBit : kBraceBit;
        body.mask = road ? kCarBit : 0;
        body.friction = 0.8f;
        body.restitution = 0;
        world.add(body);
    }
    // Hinges. At a cliff point every beam end is pinned to the world; at a free point the beams that
    // meet there are pinned to the first of them.
    for (int n = 0; n < kColumns * kRows; ++n) {
        int first = -1;
        for (size_t i = 0; i < beams.size(); ++i) {
            if (beams[i].a != n && beams[i].b != n) continue;
            const int body = static_cast<int>(first_beam_ + i);
            if (is_anchor(n)) world.add_joint(phys::Joint::revolute(world.bodies, -1, body, node_position(n)));
            else if (first < 0) first = body;
            else world.add_joint(phys::Joint::revolute(world.bodies, first, body, node_position(n)));
        }
    }
    // A structure of hinges is stiff only as far as the solver converges, so give it more sweeps.
    world.solver.iterations = 20;
    world.solver.joint_position_iterations = 8;
}

void Bridge::evaluate() {
    const Vec2 p = chassis().pos;
    if (p.y < kFallY) {
        fail();
    } else if (status == Status::Playing && p.x >= kFinishX && p.y > kCliffTop) {
        win();
        phys::scenes::set_car_throttle(world, car, 0, 40.0f);  // brake
    } else if (sim_time > 25.0) {
        fail();
    }
}

}  // namespace puzzle
