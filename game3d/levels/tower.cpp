#include "levels/tower.hpp"

#include <algorithm>
#include <cstdio>

namespace puzzle3d {

std::span<const char* const> Tower::help() const {
    static constexpr const char* kLines[] = {
        "Drag pieces onto the pad and stack them. While you hold",
        "one it hangs where it would land; let go to drop it.",
        "Q and E turn it. F stands it on end, or lays it down.",
        "Space switches physics on. The top of the tower must",
        "stay above the line for five seconds, then come to rest.",
        "The ice block is slippery; the heavy one is dense.",
    };
    return kLines;
}

const char* Tower::hint() const {
    if (status == Status::Won) return "Solved. Space to go back and build it taller.";
    if (status == Status::Failed) return "The tower did not stay above the line.";
    if (phase == Phase::Setup) return holding() ? "Q / E turn it, F stands it on end. Let go to drop it." : "Drag pieces onto the pad, then Space to test it.";
    return height() >= kGoal ? "Holding... it has to stand for five seconds." : "Below the line.";
}

std::string Tower::stats() const {
    char buf[96];
    if (phase == Phase::Running && status == Status::Playing)
        std::snprintf(buf, sizeof(buf), "Height %.1f m   Held %.1f/%.0f s", static_cast<double>(height()), std::min(held(), kHoldTime), kHoldTime);
    else
        std::snprintf(buf, sizeof(buf), "%d/%d pieces   %.1f m   best %.1f m", placed_count(), static_cast<int>(pieces.size()),
                      static_cast<double>(height()), static_cast<double>(best_height));
    return buf;
}

std::string Tower::result() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.2f m with %d piece%s", static_cast<double>(height()), placed_count(), placed_count() == 1 ? "" : "s");
    return buf;
}

void Tower::overlay(Overlay& out, const Ray* pointer) const {
    PieceLevel::overlay(out, pointer);
    // The line to beat: a dashed square at the goal height over the pad, and a post at one corner to read it by.
    const Real y = kGroundY + kGoal, h = kPad + 0.3f;
    const Vec3 c[4] = {{-h, y, -h}, {h, y, -h}, {h, y, h}, {-h, y, h}};
    for (int i = 0; i < 4; ++i) out.lines.push_back({c[i], c[(i + 1) % 4], Role::Target, true});
    out.lines.push_back({{h, kGroundY, h}, {h, y, h}, Role::Target, false});
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.1f m", static_cast<double>(kGoal));
    out.labels.push_back({{h, y + 0.35f, h}, buf, Role::Target});
    for (const Piece& p : pieces)
        if (!p.placed && (p.kind.density != 1 || p.kind.friction < 0.3f))
            out.labels.push_back({p.home + Vec3{0, p.kind.half.y + 0.4f, 0}, p.kind.label, Role::Guide});
}

Real Tower::height() const { return top_of(world, placed_bodies()) - kGroundY; }

std::vector<Piece> Tower::initial_pieces() const {
    const PieceKind block{"block", {0.6f, 0.6f, 0.6f}}, plank{"plank", {1.2f, 0.2f, 0.4f}};
    const PieceKind heavy{"heavy", {0.4f, 0.4f, 0.4f}, 0, 5.0f}, ice{"ice", {0.5f, 0.5f, 0.5f}, 0, 1.0f, 0.05f};
    std::vector<Piece> p;
    // Laid out beside the pad, on the ground, toward the camera.
    const Real x0 = -6.0f;
    for (int i = 0; i < 3; ++i) p.push_back({plank, {x0, kGroundY + 0.201f, -2.0f + 1.0f * static_cast<Real>(i)}, {}});
    p.push_back({block, {x0 + 2.6f, kGroundY + 0.601f, -1.6f}, {}});
    p.push_back({block, {x0 + 2.6f, kGroundY + 0.601f, 0.0f}, {}});
    p.push_back({heavy, {x0 + 2.6f, kGroundY + 0.401f, 1.6f}, {}});
    p.push_back({ice, {x0, kGroundY + 0.501f, 1.6f}, {}});
    return p;
}

void Tower::build_scenery() {
    add_ground(world);
    standing_ = collapsed_ = Hold{};
}

void Tower::evaluate() {
    const std::vector<int> bodies = placed_bodies();
    const bool tall = height() >= kGoal;
    // Above the line for the whole hold time, and at rest when it is judged.
    if (standing_.update(tall, kTimeStep, kHoldTime) && all_slow(world, bodies)) {
        best_height = std::max(best_height, height());
        win();
        return;
    }
    if (collapsed_.update(!tall && all_slow(world, bodies), kTimeStep, 1.0) || sim_time > 30.0) fail();
}

}  // namespace puzzle3d
