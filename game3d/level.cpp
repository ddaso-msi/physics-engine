#include "level.hpp"

namespace puzzle3d {

void add_ground(World& w) {
    Body ground = Body::fixed_box({30, 0.5f, 30}, {0, kGroundY - 0.5f, 0});
    ground.friction = 0.6f;
    ground.restitution = 0;
    w.add(ground);
}

Role Level::role(size_t body) const {
    if (body == 0) return Role::Ground;
    return world.bodies[body].type == phys3d::BodyType::Static ? Role::Scenery : Role::Tool;
}

void Level::overlay(Overlay&, const Ray*) const {}

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

}  // namespace puzzle3d
