#include "levels/bridge.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace puzzle3d {

namespace {
constexpr Real kReach = 0.5f;          // how close the pointer's ray has to pass to a connection point
constexpr Vec3 kAcross{0, 0, 1};       // every hinge's axis: across the road
// The engine has no collision filtering, so beams that meet at a point must not overlap there: the bracing
// bars stop this far short of each end (their hinges are still at the points), and the deck a hair short.
constexpr Real kBraceGap = 0.35f, kDeckGap = 0.02f;
}  // namespace

std::span<const char* const> Bridge::help() const {
    static constexpr const char* kLines[] = {
        "Drag from one connection point to a neighbour to add",
        "a beam. Click a beam to take it away. Twelve at most.",
        "Big points are fixed in rock. The right cliff's edge",
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

Vec3 Bridge::node_position(int n) {
    return {kLeftEdge + kSpacing * static_cast<Real>(n / kRows), n % kRows == 1 ? kRoadY : kLowY, 0};
}

bool Bridge::is_anchor(int n) {
    const int column = n / kRows, row = n % kRows;
    return column == 0 || (column == kColumns - 1 && row == 0);  // the right cliff's edge is crumbling
}

bool Bridge::can_connect(int a, int b) {
    if (a == b || a < 0 || b < 0 || a >= kColumns * kRows || b >= kColumns * kRows) return false;
    if (is_anchor(a) && is_anchor(b)) return false;
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

int Bridge::node_near(const Ray& ray) {
    int best = -1;
    Real best_miss = kReach;
    for (int n = 0; n < kColumns * kRows; ++n) {
        const Vec3 to = node_position(n) - ray.origin;
        const Real along = dot(to, ray.dir);
        if (along <= 0) continue;
        const Real miss = (to - ray.dir * along).length();
        if (miss < best_miss) {
            best_miss = miss;
            best = n;
        }
    }
    return best;
}

int Bridge::beam_under(const Ray& ray) const {
    const std::optional<RayHit> hit = raycast(world, ray);
    if (!hit || hit->body < static_cast<int>(first_beam_) || hit->body >= static_cast<int>(first_beam_ + beams.size())) return -1;
    return hit->body - static_cast<int>(first_beam_);
}

void Bridge::press(const Ray& ray) {
    if (phase != Phase::Setup) return;
    pointer_ = ray;
    from_ = node_near(ray);
    if (from_ >= 0) return;
    const int beam = beam_under(ray);
    if (beam >= 0) remove_beam(beams[static_cast<size_t>(beam)].a, beams[static_cast<size_t>(beam)].b);
}

void Bridge::release(const Ray& ray) {
    if (from_ < 0) return;
    const int from = from_;
    from_ = -1;
    add_beam(from, node_near(ray));
}

bool Bridge::cancel() {
    if (from_ < 0) return false;
    from_ = -1;
    return true;
}

Role Bridge::role(size_t body) const {
    if (body == 0) return Role::Ground;
    if (body >= first_beam_) return Role::Tool;
    if (world.bodies[body].type == phys3d::BodyType::Static) return Role::Scenery;
    return status == Status::Won ? Role::Done : Role::Load;
}

void Bridge::overlay(Overlay& out, const Ray* pointer) const {
    out.lines.push_back({{kFinishX, kCliffTop, -kDeckHalfWidth}, {kFinishX, kCliffTop, kDeckHalfWidth}, Role::Target, false});
    out.lines.push_back({{kFinishX, kCliffTop, kDeckHalfWidth}, {kFinishX, kCliffTop + 1.8f, kDeckHalfWidth}, Role::Target, false});
    out.labels.push_back({{kFinishX, kCliffTop + 2.1f, kDeckHalfWidth}, "finish", Role::Target});
    if (phase != Phase::Setup) return;
    for (int n = 0; n < kColumns * kRows; ++n) {
        const bool lit = n == from_ || (from_ >= 0 && can_connect(from_, n));
        out.dots.push_back({node_position(n), is_anchor(n) ? 0.2f : 0.13f, lit ? Role::Selected : is_anchor(n) ? Role::Target : Role::Guide});
    }
    if (from_ >= 0 && pointer) {
        // A line from the point to where the ray crosses the plane of the bridge.
        const Ray& r = *pointer;
        if (std::fabs(r.dir.z) > 1e-4f) {
            const Real t = -r.origin.z / r.dir.z;
            if (t > 0) out.lines.push_back({node_position(from_), r.at(t), Role::Selected, true});
        }
    }
}

void Bridge::set_throttle(Real speed, Real torque) {
    for (int d : car.drives) {
        phys3d::Joint& j = world.joints[static_cast<size_t>(d)];
        j.motor_speed = -speed;  // a wheel rolling toward +x turns clockwise about +z, i.e. negatively
        j.max_motor = torque;
    }
}

void Bridge::build() {
    add_ground(world);
    const Real right_edge = kLeftEdge + kSpacing * static_cast<Real>(kColumns - 1);
    const Real h = (kCliffTop - kGroundY) * 0.5f;
    world.add(Body::fixed_box({6, h, 4}, {kLeftEdge - 6, kGroundY + h, 0}));
    world.add(Body::fixed_box({6, h, 4}, {right_edge + 6, kGroundY + h, 0}));

    // The car: a body on four wheels, each on a motorised hinge.
    constexpr Real kWheel = 0.42f;
    const Vec3 at{kLeftEdge - 3.0f, kCliffTop + kWheel + 0.38f + 0.01f, 0};
    Body body = Body::solid_box({1.1f, 0.22f, 0.75f}, 1.5f, at);
    body.friction = 0.4f;
    car.chassis = world.add(body);
    int k = 0;
    for (const Real dx : {-0.8f, 0.8f})
        for (const Real dz : {-0.9f, 0.9f}) {
            Body wheel = Body::solid_sphere(kWheel, 1.5f, at + Vec3{dx, -0.38f, dz});
            wheel.friction = 1.2f;
            wheel.restitution = 0.05f;
            car.wheels[k] = world.add(wheel);
            phys3d::Joint axle = phys3d::Joint::hinge(world.bodies, car.chassis, car.wheels[k], wheel.pos, kAcross);
            axle.enable_motor = true;
            car.drives[k] = world.add_joint(axle);
            ++k;
        }
    set_throttle(7.0f, 20.0f);

    // Beams: the road as a deck the width of the road, the bracing as one bar under the middle.
    first_beam_ = world.bodies.size();
    for (const Beam& beam : beams) {
        const Vec3 a = node_position(beam.a), b = node_position(beam.b), d = b - a;
        const bool road = is_road(beam);
        const Real half_len = d.length() * 0.5f - (road ? kDeckGap : kBraceGap);
        const Quat turn = Quat::from_axis_angle(kAcross, std::atan2(d.y, d.x));
        Body body_ = Body::solid_box({half_len, road ? 0.1f : 0.07f, road ? kDeckHalfWidth : 0.07f}, road ? 2.0f : 6.0f, (a + b) * 0.5f, turn);
        body_.friction = 0.8f;
        body_.restitution = 0;
        world.add(body_);
    }
    // Hinges. At a fixed point every beam end is pinned to the world; at a free point the beams that meet
    // there are pinned to one of them, a road beam if there is one, so neighbouring deck planks are joined
    // to each other and do not rub.
    for (int n = 0; n < kColumns * kRows; ++n) {
        std::vector<int> here;
        for (size_t i = 0; i < beams.size(); ++i)
            if (beams[i].a == n || beams[i].b == n) here.push_back(static_cast<int>(i));
        std::stable_partition(here.begin(), here.end(), [&](int i) { return is_road(beams[static_cast<size_t>(i)]); });
        for (size_t i = 0; i < here.size(); ++i) {
            const int body_index = static_cast<int>(first_beam_) + here[i];
            if (is_anchor(n)) world.add_joint(phys3d::Joint::hinge(world.bodies, -1, body_index, node_position(n), kAcross));
            else if (i > 0) world.add_joint(phys3d::Joint::hinge(world.bodies, static_cast<int>(first_beam_) + here[0], body_index, node_position(n), kAcross));
        }
    }
    world.solver.iterations = 20;
}

void Bridge::evaluate() {
    const Vec3 p = chassis().pos;
    if (p.y < kFallY) {
        fail();
    } else if (status == Status::Playing && p.x >= kFinishX && p.y > kCliffTop) {
        win();
        set_throttle(0, 40.0f);  // brake
    } else if (sim_time > 25.0) {
        fail();
    }
}

}  // namespace puzzle3d
