#include "levels/chain_reaction.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace puzzle3d {

namespace {
// Angle round the curve, from +x (0) to +z (pi/2), of slot k.
Real slot_angle(int k) { return static_cast<Real>(k) * ChainReaction::kPitch / ChainReaction::kRadius + 0.05f; }
}  // namespace

Vec3 ChainReaction::slot_position(int k) {
    const Real a = slot_angle(k);
    return {kRadius * std::cos(a), kGroundY + kDominoHalf.y + 0.001f, kRadius * std::sin(a)};
}

// A domino's own x axis is its thin side, so it should point along the curve: (-sin a, 0, cos a). A turn by
// yaw about +y takes x to (cos yaw, 0, -sin yaw), so yaw = -(a + pi/2).
Real ChainReaction::slot_yaw(int k) { return -(slot_angle(k) + 3.14159265f / 2); }

std::span<const char* const> ChainReaction::help() const {
    static constexpr const char* kLines[] = {
        "Drag the spare dominoes into the marked gap. While you",
        "hold one it hangs where it would land; let go to drop it.",
        "Q and E turn it: a domino falls the way it faces, so",
        "turn each one to follow the curve.",
        "Space starts the ball rolling. The chain has to reach",
        "the far ball and tip it into the bin.",
    };
    return kLines;
}

const char* ChainReaction::hint() const {
    if (status == Status::Won) return "Solved.";
    if (status == Status::Failed) return "The chain stopped before the ball reached the bin.";
    if (phase == Phase::Setup) return holding() ? "Q / E turn it to follow the curve. Let go to drop it." : "Drag dominoes into the gap, then Space to start.";
    return "Space to stop and rearrange.";
}

std::string ChainReaction::stats() const {
    char buf[48];
    std::snprintf(buf, sizeof(buf), "Dominoes placed %d/%d", placed_count(), static_cast<int>(pieces.size()));
    return buf;
}

std::string ChainReaction::result() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%d spare domino%s, %.1f s", placed_count(), placed_count() == 1 ? "" : "es", win_time);
    return buf;
}

void ChainReaction::overlay(Overlay& out, const Ray* pointer) const {
    out.zones.push_back({kZoneLo, kZoneHi, Role::Target});
    out.labels.push_back({{(kZoneLo.x + kZoneHi.x) * 0.5f, kZoneHi.y + 0.5f, (kZoneLo.z + kZoneHi.z) * 0.5f}, "bin", Role::Target});
    // The curve the dominoes follow, dashed on the ground, so the gap reads as part of it.
    for (int k = 0; k + 1 < kSlots; ++k) {
        Vec3 a = slot_position(k), b = slot_position(k + 1);
        a.y = b.y = kGroundY + 0.02f;
        out.lines.push_back({a, b, Role::Guide, true});
    }
    PieceLevel::overlay(out, pointer);
}

std::vector<Piece> ChainReaction::initial_pieces() const {
    std::vector<Piece> p;
    for (int i = 0; i < kGapCount; ++i) {
        Piece d;
        d.kind = {"domino", kDominoHalf, 0, 2.0f, 0.6f, 0.05f};
        d.home = {1.0f + 0.8f * static_cast<Real>(i), kGroundY + kDominoHalf.y + 0.001f, 1.2f};
        p.push_back(d);
    }
    return p;
}

phys3d::AABB ChainReaction::build_area() const {
    // A box around the empty slots, no higher than a standing domino: so a spare can stand anywhere in the
    // gap, turned any way, but not on top of anything.
    phys3d::AABB box{slot_position(kGapFirst), slot_position(kGapFirst)};
    for (int k = kGapFirst; k < kGapFirst + kGapCount; ++k) {
        const Vec3 p = slot_position(k);
        box.lo = {std::min(box.lo.x, p.x), kGroundY, std::min(box.lo.z, p.z)};
        box.hi = {std::max(box.hi.x, p.x), kGroundY + 2 * kDominoHalf.y + 0.05f, std::max(box.hi.z, p.z)};
    }
    box.lo = box.lo - Vec3{0.5f, 0, 0.5f};
    box.hi = box.hi + Vec3{0.5f, 0, 0.5f};
    return box;
}

Role ChainReaction::scenery_role(size_t body) const {
    if (static_cast<int>(body) == starter) return Role::Load;
    if (static_cast<int>(body) == target) return phase == Phase::Running && target_in_bin() ? Role::Done : Role::Target;
    if (world.bodies[body].type == phys3d::BodyType::Dynamic) return Role::Guide;
    return PieceLevel::scenery_role(body);
}

void ChainReaction::build_scenery() {
    add_ground(world);
    chain_.clear();
    // The ramp: it slopes down toward +z and ends at the ground just before the first domino.
    const Real slope = 0.42f, half_len = 2.2f;
    const Quat tilt = Quat::from_axis_angle({1, 0, 0}, slope);  // turning about +x lowers the +z end
    const Vec3 along = phys3d::rotate(tilt, Vec3{0, 0, 1}), up = phys3d::rotate(tilt, Vec3{0, 1, 0});
    const Vec3 first = slot_position(0);
    const Vec3 low_end{first.x, kGroundY, first.z - 0.9f};
    const Vec3 centre = low_end - along * half_len + up * -0.1f;
    world.add(Body::fixed_box({0.6f, 0.1f, half_len}, centre, tilt));
    Body ball = Body::solid_sphere(0.3f, 3.0f, centre - along * (half_len - 0.5f) + up * 0.402f);
    ball.friction = 0.6f;
    ball.restitution = 0.05f;
    starter = world.add(ball);

    for (int k = 0; k < kSlots; ++k) {
        if (k >= kGapFirst && k < kGapFirst + kGapCount) continue;
        Body d = Body::solid_box(kDominoHalf, 2.0f, slot_position(k), Quat::from_axis_angle({0, 1, 0}, slot_yaw(k)));
        d.friction = 0.6f;
        d.restitution = 0.05f;
        chain_.push_back(world.add(d));
    }

    // The last ball, just past the last domino, and the bin beyond it. Slippery, so the last domino cannot
    // wedge it against the ground.
    const Real a_end = slot_angle(kSlots - 1) + 0.45f / kRadius;
    Body goal = Body::solid_sphere(0.3f, 1.0f, {kRadius * std::cos(a_end), kGroundY + 0.301f, kRadius * std::sin(a_end)});
    goal.friction = 0.1f;
    goal.restitution = 0.1f;
    target = world.add(goal);
    chain_.push_back(target);
    const Vec3 mid = (kZoneLo + kZoneHi) * 0.5f;
    const Real half_x = (kZoneHi.x - kZoneLo.x) * 0.5f, half_z = (kZoneHi.z - kZoneLo.z) * 0.5f;
    // The ball comes in from +z, so the bin is open on that side: a wall at the back and one down each side.
    world.add(Body::fixed_box({half_x + 0.2f, 0.25f, 0.1f}, {mid.x, 0.25f, kZoneLo.z - 0.1f}));
    world.add(Body::fixed_box({0.1f, 0.25f, half_z}, {kZoneLo.x - 0.1f, 0.25f, mid.z}));
    world.add(Body::fixed_box({0.1f, 0.25f, half_z}, {kZoneHi.x + 0.1f, 0.25f, mid.z}));
    resting_ = stalled_ = Hold{};
}

void ChainReaction::evaluate() {
    if (resting_.update(target_in_bin(), kTimeStep, kRestTime)) {
        win();
        return;
    }
    // Lost once the chain has stopped: every domino and the last ball at rest (the first ball may roll on).
    std::vector<int> moving = chain_;
    for (const Piece& p : pieces)
        if (p.body >= 0) moving.push_back(p.body);
    if ((stalled_.update(all_slow(world, moving), kTimeStep, 1.0) && sim_time > 3.0) || sim_time > 25.0) fail();
}

}  // namespace puzzle3d
