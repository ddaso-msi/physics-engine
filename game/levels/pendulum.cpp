#include "levels/pendulum.hpp"

#include <scenes.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace puzzle {

namespace {
Vec2 bob_position(Real angle) {
    return Pendulum::kPivot + Vec2{Pendulum::kRod * std::sin(angle), -Pendulum::kRod * std::cos(angle)};
}
constexpr Real kDegree = phys::kPi / 180.0f;
}  // namespace

std::span<const char* const> Pendulum::help() const {
    static constexpr const char* kLines[] = {
        "Drag the weight along its arc to pull it back.",
        "Q and E change the angle one degree at a time.",
        "Space lets go. The crate has to land in the bin",
        "and stay there. Space again to try another angle.",
    };
    return kLines;
}

const char* Pendulum::hint() const {
    if (status == Status::Won) return "Solved.";
    if (status == Status::Failed) return crate_in_bin() ? "The crate did not stay in the bin." : "The crate missed the bin.";
    if (phase == Phase::Setup) return "Drag the weight to set the angle, then Space to release.";
    return "Space to stop and change the angle.";
}

std::string Pendulum::stats() const {
    char buf[48];
    std::snprintf(buf, sizeof(buf), "Angle %.0f deg", static_cast<double>(-angle / kDegree));
    return buf;
}

std::string Pendulum::result() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.0f degrees, %.1f s", static_cast<double>(-angle / kDegree), win_time);
    return buf;
}

Role Pendulum::role(size_t body) const {
    if (static_cast<int>(body) == bob) return dragging_ ? Role::Selected : Role::Tool;
    if (static_cast<int>(body) == crate) return phase == Phase::Running && crate_in_bin() ? Role::Done : Role::Target;
    return Role::Scenery;
}

void Pendulum::overlay(Overlay& out, Vec2, bool) const {
    out.zones.push_back({kZoneLo, kZoneHi, Role::Target});
    out.labels.push_back({{kZoneLo.x + 0.45f, kZoneHi.y + 0.5f}, "bin", Role::Target});
    if (phase != Phase::Setup) return;
    for (Real a = kMinAngle; a <= kMaxAngle + 0.001f; a += 0.075f) out.dots.push_back({bob_position(a), 0.04f, Role::Guide});
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.0f deg", static_cast<double>(-angle / kDegree));
    out.labels.push_back({bob_position(angle) + Vec2{-0.6f, 1.1f}, buf, Role::Tool});
}

void Pendulum::set_angle(Real a) {
    if (phase != Phase::Setup) return;
    angle = std::clamp(a, kMinAngle, kMaxAngle);
    rebuild();
}

void Pendulum::press(Vec2 p) {
    if (phase != Phase::Setup || distance(p, world.bodies[static_cast<size_t>(bob)].pos) > kBobRadius + 0.5f) return;
    dragging_ = true;
    angle_before_drag_ = angle;
}

void Pendulum::drag(Vec2 p) {
    if (!dragging_) return;
    set_angle(std::atan2(p.x - kPivot.x, kPivot.y - p.y));
}

void Pendulum::rotate(int direction) { set_angle(angle - static_cast<Real>(direction) * kDegree); }

bool Pendulum::cancel() {
    if (!dragging_) return false;
    dragging_ = false;
    set_angle(angle_before_drag_);
    return true;
}

void Pendulum::build() {
    using phys::scenes::add_static_box;
    add_room(world);
    add_static_box(world, {1.6f, 0.15f}, {8.5f, 3.65f});                    // the table: top at 3.8
    add_static_box(world, {0.15f, (3.5f - kFloorTop) * 0.5f}, {8.5f, kFloorTop + (3.5f - kFloorTop) * 0.5f});
    add_static_box(world, {0.1f, 0.3f}, {kZoneLo.x - 0.1f, kFloorTop + 0.3f});  // the bin's two lips
    add_static_box(world, {0.1f, 0.3f}, {kZoneHi.x + 0.1f, kFloorTop + 0.3f});

    Body box(phys::Shape::make_polygon(phys::Polygon::box(0.4f, 0.4f)), {7.55f, 3.8f + 0.401f}, 0);
    box.friction = 0.5f;
    box.restitution = 0.1f;
    crate = world.add(box);

    // The rod is a rigid distance joint to the world, made with the weight already at the chosen angle.
    Body weight(phys::Shape::make_circle(kBobRadius), bob_position(angle), 0, phys::BodyType::Dynamic, 4.0f);
    weight.friction = 0.4f;
    weight.restitution = 0.2f;
    weight.allow_sleep = false;  // it may hang still at the top of a swing; it must not doze off there
    bob = world.add(weight);
    world.add_joint(phys::Joint::distance(world.bodies, -1, bob, kPivot, weight.pos));
    resting_ = stalled_ = Hold{};
}

void Pendulum::evaluate() {
    const int the_crate[] = {crate};
    const bool slow = all_slow(world, the_crate);
    const bool in_bin = crate_in_bin();
    if (resting_.update(in_bin && slow, kTimeStep, kRestTime)) {
        win();
        return;
    }
    // Lost once the crate has come to rest somewhere else, with time allowed for the weight to arrive.
    const bool stuck = stalled_.update(slow && !in_bin, kTimeStep, 1.5) && sim_time > 4.0;
    if (stuck || sim_time > 20.0) fail();
}

}  // namespace puzzle
