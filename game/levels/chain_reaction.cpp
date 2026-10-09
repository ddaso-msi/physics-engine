#include "levels/chain_reaction.hpp"

#include <scenes.hpp>

#include <cmath>
#include <cstdio>

namespace puzzle {

std::span<const char* const> ChainReaction::help() const {
    static constexpr const char* kLines[] = {
        "Drag the spare dominoes from the tray into the",
        "marked gap. They drop onto the shelf where you let go.",
        "Q and E tilt the one you are holding.",
        "Space starts the ball rolling. The chain has to carry",
        "all the way to the far ball and tip it into the bin.",
    };
    return kLines;
}

const char* ChainReaction::hint() const {
    if (status == Status::Won) return "Solved.";
    if (status == Status::Failed) return "The chain stopped before the ball reached the bin.";
    if (phase == Phase::Setup) return holding() ? "Let go inside the marked gap." : "Drag dominoes into the gap, then Space to start.";
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

void ChainReaction::overlay(Overlay& out, Vec2 pointer, bool pointer_active) const {
    out.zones.push_back({kZoneLo, kZoneHi, Role::Target});
    out.labels.push_back({{kZoneLo.x + 0.9f, kZoneHi.y + 0.5f}, "bin", Role::Target});
    PieceLevel::overlay(out, pointer, pointer_active);
}

std::vector<Piece> ChainReaction::initial_pieces() const {
    std::vector<Piece> p;
    for (int i = 0; i < 3; ++i) {
        Piece d;
        d.kind = {"domino", kDominoHalf, 0, 2.0f, 0.6f, 0.05f};
        d.home = {7.3f + 0.7f * static_cast<Real>(i), 6.6f};
        p.push_back(d);
    }
    return p;
}

Role ChainReaction::scenery_role(size_t body) const {
    if (static_cast<int>(body) == starter) return Role::Load;
    if (static_cast<int>(body) == target) return phase == Phase::Running && target_in_bin() ? Role::Done : Role::Target;
    return PieceLevel::scenery_role(body);
}

void ChainReaction::build_scenery() {
    using phys::scenes::add_static_box;
    add_room(world);
    const Real mid = (kShelfLeft + kShelfRight) * 0.5f, half = (kShelfRight - kShelfLeft) * 0.5f;
    add_static_box(world, {half, 0.15f}, {mid, kShelfTop - 0.15f});
    for (const Real x : {kShelfLeft + 1.0f, kShelfRight - 1.0f}) {
        const Real h = (kShelfTop - 0.3f - kFloorTop) * 0.5f;
        add_static_box(world, {0.15f, h}, {x, kFloorTop + h});
    }
    add_static_box(world, {0.1f, 0.3f}, {kZoneLo.x - 0.1f, kFloorTop + 0.3f});  // the bin's lips
    add_static_box(world, {0.1f, 0.3f}, {kZoneHi.x + 0.1f, kFloorTop + 0.3f});

    // The ramp slopes down to the right and meets the shelf's left end.
    constexpr Real kSlope = -0.5f, kRampHalf = 1.7f;
    const Vec2 along{std::cos(kSlope), std::sin(kSlope)}, up{-along.y, along.x};
    const Vec2 ramp_centre = Vec2{kShelfLeft + 0.4f, kShelfTop + 0.15f} - along * kRampHalf;
    add_static_box(world, {kRampHalf, 0.1f}, ramp_centre, kSlope);

    Body ball(phys::Shape::make_circle(0.3f), ramp_centre - along * 1.2f + up * 0.402f, 0, phys::BodyType::Dynamic, 3.0f);
    ball.friction = 0.6f;
    ball.restitution = 0.05f;
    starter = world.add(ball);

    auto domino = [&](Real x) {
        Body d(phys::Shape::make_polygon(phys::Polygon::box(kDominoHalf.x, kDominoHalf.y)), {x, kShelfTop + kDominoHalf.y + 0.001f}, 0,
               phys::BodyType::Dynamic, 2.0f);
        d.friction = 0.6f;
        d.restitution = 0.05f;
        world.add(d);
    };
    for (const Real x : {kGapLeft - 1.3f, kGapLeft - 0.65f, kGapLeft}) domino(x);
    for (int i = 0; i < 5; ++i) domino(kGapRight + 0.65f * static_cast<Real>(i));

    Body goal(phys::Shape::make_circle(0.3f), {kGapRight + 0.65f * 4 + 0.42f, kShelfTop + 0.301f}, 0);
    // Slippery on purpose. With ordinary friction the last domino leans on the ball and wedges it against
    // the shelf (a doorstop), and the chain ends one step short.
    goal.friction = 0.1f;
    goal.restitution = 0.1f;
    target = world.add(goal);
    resting_ = stalled_ = Hold{};
}

void ChainReaction::evaluate() {
    if (resting_.update(target_in_bin(), kTimeStep, kRestTime)) {
        win();
        return;
    }
    // Lost when nothing is moving any more (after giving the first ball time to get going).
    const bool quiet = world.stats().awake_bodies == 0;
    if ((stalled_.update(quiet, kTimeStep, 0.5) && sim_time > 2.0) || sim_time > 25.0) fail();
}

}  // namespace puzzle
