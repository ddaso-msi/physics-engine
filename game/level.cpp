#include "level.hpp"

#include <scenes.hpp>

namespace puzzle {

void add_room(World& w) {
    using phys::scenes::add_static_box;
    const Real t = kFloorTop * 0.5f;
    add_static_box(w, {kWorldW * 0.5f, t}, {kWorldW * 0.5f, t});
    add_static_box(w, {t, kWorldH * 0.5f}, {t, kWorldH * 0.5f});
    add_static_box(w, {t, kWorldH * 0.5f}, {kWorldW - t, kWorldH * 0.5f});
    add_static_box(w, {kWorldW * 0.5f, t}, {kWorldW * 0.5f, kWorldH + t});
}

Role Level::role(size_t body) const {
    return world.bodies[body].type == phys::BodyType::Static ? Role::Scenery : Role::Tool;
}

void Level::overlay(Overlay&, Vec2, bool) const {}

void Level::rebuild() {
    world = World{};
    build();
}

void Level::begin_run() {
    sim_time = win_time = 0;
    status = Status::Playing;
    if (tried_) ++attempt;
    tried_ = false;
}

void Level::reset() {
    cancel();
    clear_setup();
    begin_run();
    phase = has_setup() ? Phase::Setup : Phase::Running;
    rebuild();
}

void Level::start() {
    if (phase != Phase::Setup) return;
    cancel();
    begin_run();
    rebuild();
    phase = Phase::Running;
    mark_tried();
}

void Level::stop() {
    if (phase != Phase::Running || !has_setup()) return;
    begin_run();
    phase = Phase::Setup;
    rebuild();
}

void Level::step(Real dt) {
    if (phase != Phase::Running) return;
    world.step(dt);
    sim_time += static_cast<double>(dt);
    if (status != Status::Won) evaluate();
}

void Level::win() {
    if (status == Status::Won) return;
    status = Status::Won;
    win_time = sim_time;
}

void Level::fail() {
    if (status == Status::Playing) status = Status::Failed;
}

}  // namespace puzzle
