#include "levels/pendulum.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace puzzle3d {

namespace {
constexpr Real kDegree = 3.14159265f / 180.0f;
}

Vec3 Pendulum::bob_position(Real t, Real h) {
    // Straight back is -x; a positive heading swings the pull round toward +z.
    const Vec3 out{-std::cos(h), 0, std::sin(h)};
    return kPivot + out * (kRod * std::sin(t)) + Vec3{0, -kRod * std::cos(t), 0};
}

std::span<const char* const> Pendulum::help() const {
    static constexpr const char* kLines[] = {
        "Drag the weight to pull it back. Further back is faster;",
        "pulling it to one side sends the ball off to the other.",
        "Q and E turn the pull one degree at a time.",
        "Space lets go. The ball has to land in the bin and",
        "stay there. Space again to try another pull.",
    };
    return kLines;
}

const char* Pendulum::hint() const {
    if (status == Status::Won) return "Solved.";
    if (status == Status::Failed) return target_in_bin() ? "The ball did not stay in the bin." : "The ball missed the bin.";
    if (phase == Phase::Setup) return "Drag the weight to pull it back, then Space to let go.";
    return "Space to stop and change the pull.";
}

std::string Pendulum::stats() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "Pull %.0f deg   Aim %+.0f deg", static_cast<double>(tilt / kDegree), static_cast<double>(heading / kDegree));
    return buf;
}

std::string Pendulum::result() const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.0f deg back, %+.0f deg aim, %.1f s", static_cast<double>(tilt / kDegree),
                  static_cast<double>(heading / kDegree), win_time);
    return buf;
}

Role Pendulum::role(size_t body) const {
    if (static_cast<int>(body) == bob) return dragging_ ? Role::Selected : Role::Tool;
    if (static_cast<int>(body) == target) return phase == Phase::Running && target_in_bin() ? Role::Done : Role::Target;
    return Level::role(body);
}

void Pendulum::overlay(Overlay& out, const Ray*) const {
    out.zones.push_back({kZoneLo, kZoneHi, Role::Target});
    out.labels.push_back({{(kZoneLo.x + kZoneHi.x) * 0.5f, kZoneHi.y + 0.6f, 0}, "bin", Role::Target});
    if (phase != Phase::Setup) return;
    // The arc the ball can be pulled along at this heading, and a dashed line on the ground showing the
    // direction the ball will swing through.
    for (Real t = kMinTilt; t <= kMaxTilt + 0.001f; t += 0.07f) out.dots.push_back({bob_position(t, heading), 0.04f, Role::Guide});
    const Vec3 dir{std::cos(heading), 0, -std::sin(heading)};
    const Vec3 under{kPivot.x, kGroundY + 0.02f, kPivot.z};
    out.lines.push_back({under - dir * 2.0f, under + dir * 9.0f, Role::Guide, true});
    char buf[48];
    std::snprintf(buf, sizeof(buf), "%.0f deg", static_cast<double>(tilt / kDegree));
    out.labels.push_back({bob_position(tilt, heading) + Vec3{0, 1.0f, 0}, buf, Role::Tool});
}

void Pendulum::set_pull(Real new_tilt, Real new_heading) {
    if (phase != Phase::Setup) return;
    tilt = std::clamp(new_tilt, kMinTilt, kMaxTilt);
    heading = std::clamp(new_heading, -0.9f, 0.9f);
    rebuild();
}

void Pendulum::press(const Ray& ray) {
    if (phase != Phase::Setup) return;
    const std::optional<RayHit> hit = raycast(world.bodies[static_cast<size_t>(bob)], ray);
    // A generous grab: the ball itself, or anywhere within half a metre of it on screen.
    const Vec3 to = world.bodies[static_cast<size_t>(bob)].pos - ray.origin;
    const Real miss = (to - ray.dir * dot(to, ray.dir)).length();
    if (!hit && miss > kBobRadius + 0.5f) return;
    dragging_ = true;
    tilt_before_ = tilt;
    heading_before_ = heading;
}

void Pendulum::drag(const Ray& ray) {
    if (!dragging_) return;
    // Where the ray meets the sphere the rod sweeps; failing that, the point of the ray closest to it. The
    // ray goes in through the near side of the sphere and out through the far side: of the two, take the one
    // nearer the weight as it is now, so the weight follows the pointer over the back of the sphere too.
    const Vec3 m = ray.origin - kPivot;
    const Real b = dot(m, ray.dir), c = m.length_sq() - kRod * kRod, disc = b * b - c;
    Vec3 p = ray.at(-b);
    if (disc >= 0) {
        const Vec3 near_side = ray.at(-b - std::sqrt(disc)), far_side = ray.at(-b + std::sqrt(disc));
        const Vec3 now = bob_position(tilt, heading);
        p = distance(near_side, now) <= distance(far_side, now) ? near_side : far_side;
    }
    const Vec3 d = p - kPivot;
    const Real horizontal = std::sqrt(d.x * d.x + d.z * d.z);
    set_pull(std::atan2(horizontal, -d.y), horizontal > 1e-4f ? std::atan2(d.z, -d.x) : heading);
}

void Pendulum::rotate(int direction) { set_pull(tilt, heading + static_cast<Real>(direction) * kDegree); }

bool Pendulum::cancel() {
    if (!dragging_) return false;
    dragging_ = false;
    set_pull(tilt_before_, heading_before_);
    return true;
}

void Pendulum::build() {
    add_ground(world);
    // The table, top at kTableTop, and its leg. The ball waits on it where the weight, swinging up the far
    // side, meets it square on; so the weight is already rising when it hits, and the aim matters.
    world.add(Body::fixed_box({0.9f, 0.15f, 0.9f}, {3.5f, kTableTop - 0.15f, 0}));
    world.add(Body::fixed_box({0.15f, (kTableTop - 0.3f) * 0.5f, 0.15f}, {3.5f, (kTableTop - 0.3f) * 0.5f, 0}));
    world.add(Body::fixed_box({0.25f, 4.05f, 0.25f}, {kPivot.x, 4.05f, -2.4f}));    // the gallows: a post beside the swing...
    world.add(Body::fixed_box({0.2f, 0.2f, 1.25f}, {kPivot.x, kPivot.y + 0.4f, -1.25f}));  // ...and an arm out over the pivot
    // The bin: four walls half a metre high, enough to stop a ball rolling out.
    const Real wall = 0.25f, mid_x = (kZoneLo.x + kZoneHi.x) * 0.5f, half_x = (kZoneHi.x - kZoneLo.x) * 0.5f + 0.2f;
    world.add(Body::fixed_box({0.1f, wall, kZoneHi.z + 0.2f}, {kZoneLo.x - 0.1f, wall, 0}));
    world.add(Body::fixed_box({0.1f, wall, kZoneHi.z + 0.2f}, {kZoneHi.x + 0.1f, wall, 0}));
    world.add(Body::fixed_box({half_x, wall, 0.1f}, {mid_x, wall, kZoneLo.z - 0.1f}));
    world.add(Body::fixed_box({half_x, wall, 0.1f}, {mid_x, wall, kZoneHi.z + 0.1f}));

    Body ball_target = Body::solid_sphere(kTargetRadius, 1, kTargetStart);
    ball_target.friction = 0.5f;
    ball_target.restitution = 0.1f;
    target = world.add(ball_target);

    Body ball = Body::solid_sphere(kBobRadius, 4.0f, bob_position(tilt, heading));
    ball.friction = 0.4f;
    ball.restitution = 0.2f;
    ball.allow_sleep = false;  // it can hang still at the top of a swing; it must not doze off there
    bob = world.add(ball);
    world.add_joint(phys3d::Joint::distance(world.bodies, -1, bob, kPivot, ball.pos));
    resting_ = stalled_ = Hold{};
}

void Pendulum::evaluate() {
    const int the_target[] = {target};
    const bool in_bin = target_in_bin();
    // Nothing in the engine slows a rolling ball, so the ball may never be quite still: in the bin, held in
    // by its walls, for long enough is the win.
    if (resting_.update(in_bin, kTimeStep, kRestTime)) {
        win();
        return;
    }
    // Lost once the ball is down on the ground somewhere else, or has stopped on the table.
    const Body& ball = world.bodies[static_cast<size_t>(target)];
    const bool landed = ball.pos.y < kGroundY + kTargetRadius + 0.05f && !in_bin;
    const bool stuck = all_slow(world, the_target) && !in_bin && sim_time > 4.0;
    if (stalled_.update(landed || stuck, kTimeStep, 1.0) || sim_time > 20.0) fail();
}

}  // namespace puzzle3d
