#include "levels/tower.hpp"

#include <algorithm>
#include <cstdio>

namespace puzzle {

std::span<const char* const> Tower::help() const {
    static constexpr const char* kLines[] = {
        "Drag pieces from the tray and stack them. A piece",
        "drops straight down onto whatever is below it.",
        "Q and E turn the piece you are holding on its side.",
        "Space switches physics on. The top of the tower must",
        "stay above the line for five seconds, then come to rest.",
        "The ice block is slippery; the heavy one is dense.",
    };
    return kLines;
}

const char* Tower::hint() const {
    if (status == Status::Won) return "Solved. Space to go back and build it taller.";
    if (status == Status::Failed) return "The tower did not stay above the line.";
    if (phase == Phase::Setup) return holding() ? "Q / E turn it. Let go to drop it." : "Drag pieces to build, then Space to test it.";
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

void Tower::overlay(Overlay& out, Vec2 pointer, bool pointer_active) const {
    PieceLevel::overlay(out, pointer, pointer_active);
    const Real y = kFloorTop + kGoal;
    out.lines.push_back({{5.6f, y}, {kWorldW - 0.6f, y}, Role::Target, true});
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.1f m", static_cast<double>(kGoal));
    out.labels.push_back({{kWorldW - 2.2f, y + 0.4f}, buf, Role::Target});
    // The two odd pieces look like plain blocks, so name them while they wait in the tray.
    for (const Piece& p : pieces)
        if (!p.placed && (p.kind.density != 1 || p.kind.friction < 0.3f))
            out.labels.push_back({p.home + Vec2{-p.kind.half.x, -p.kind.half.y - 0.1f}, p.kind.label, Role::Guide});
}

Real Tower::height() const { return top_of(world, placed_bodies()) - kFloorTop; }

std::vector<Piece> Tower::initial_pieces() const {
    const PieceKind block{"block", {0.6f, 0.6f}}, plank{"plank", {1.2f, 0.2f}};
    const PieceKind heavy{"heavy", {0.4f, 0.4f}, 0, 5.0f}, ice{"ice", {0.5f, 0.5f}, 0, 1.0f, 0.05f};
    std::vector<Piece> p;
    for (const Real y : {8.9f, 8.3f, 7.7f}) p.push_back({plank, {2.9f, y}, {}});
    p.push_back({block, {2.1f, 6.5f}, {}});
    p.push_back({block, {3.7f, 6.5f}, {}});
    p.push_back({heavy, {2.1f, 4.9f}, {}});
    p.push_back({ice, {3.8f, 5.0f}, {}});
    return p;
}

void Tower::build_scenery() {
    add_room(world);
    standing_ = collapsed_ = Hold{};
}

void Tower::evaluate() {
    const std::vector<int> bodies = placed_bodies();
    const bool tall = height() >= kGoal;
    // Above the line for the whole hold time, and at rest when it is judged: a tower that is still
    // swaying after five seconds has to finish swaying first.
    if (standing_.update(tall, kTimeStep, kHoldTime) && all_slow(world, bodies)) {
        best_height = std::max(best_height, height());
        win();
        return;
    }
    if ((collapsed_.update(!tall && all_slow(world, bodies), kTimeStep, 1.0)) || sim_time > 30.0) fail();
}

}  // namespace puzzle
