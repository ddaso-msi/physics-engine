#include "scenes.hpp"

#include <cmath>

namespace phys::scenes {

namespace {

Body make_box(Vec2 half, Vec2 pos, Real angle, Real density = 1) {
    return Body(Shape::make_polygon(Polygon::box(half.x, half.y)), pos, angle, BodyType::Dynamic, density);
}

// All parts of one assembly share a category bit and refuse to collide with it.
void set_group(Body& b, int group) {
    b.category = static_cast<std::uint16_t>(1u << group);
    b.mask = static_cast<std::uint16_t>(0xFFFFu & ~(1u << group));
}

Joint& add_limited_hinge(World& w, int parent, int child, Vec2 anchor, Real lower, Real upper, std::vector<int>* record) {
    Joint j = Joint::revolute(w.bodies, parent, child, anchor);
    j.enable_limit = true;
    j.lower = lower;
    j.upper = upper;
    const int index = w.add_joint(j);
    if (record) record->push_back(index);
    return w.joints[static_cast<size_t>(index)];
}

}  // namespace

int add_static_box(World& w, Vec2 half_extents, Vec2 pos, Real angle) {
    Body b(Shape::make_polygon(Polygon::box(half_extents.x, half_extents.y)), pos, angle, BodyType::Static);
    b.restitution = 0;
    b.friction = 0.6f;
    return w.add(b);
}

Cradle build_cradle(World& w, Vec2 top_center, int count, Real rope, Real r) {
    Cradle c;
    c.rope = rope;
    for (int i = 0; i < count; ++i) {
        // A hair more than a diameter apart so the balls do not start out overlapping.
        const Real x = top_center.x + (static_cast<Real>(i) - static_cast<Real>(count - 1) * 0.5f) * (2 * r + 0.001f);
        Body ball(Shape::make_circle(r), {x, top_center.y - rope}, 0, BodyType::Dynamic, 3);
        ball.restitution = 1;
        ball.friction = 0;
        const int idx = w.add(ball);
        c.balls.push_back(idx);
        c.pivots.push_back({x, top_center.y});
        w.add_joint(Joint::distance(w.bodies, -1, idx, {x, top_center.y}, ball.pos));
    }
    return c;
}

void lift_cradle_ball(World& w, const Cradle& c, int which, Real angle) {
    Body& ball = w.bodies[static_cast<size_t>(c.balls[static_cast<size_t>(which)])];
    const Vec2 pivot = c.pivots[static_cast<size_t>(which)];
    ball.pos = pivot + Vec2{c.rope * std::sin(angle), -c.rope * std::cos(angle)};
    ball.vel = {};
    ball.wake();
}

Ragdoll build_ragdoll(World& w, Vec2 torso, int group) {
    Ragdoll d;
    auto part = [&](Body b) {
        set_group(b, group);
        b.friction = 0.6f;
        b.restitution = 0.05f;
        return w.add(b);
    };

    d.torso = part(make_box({0.28f, 0.55f}, torso, 0));
    const Vec2 neck = torso + Vec2{0, 0.55f};
    d.head = part(Body(Shape::make_circle(0.28f), neck + Vec2{0, 0.25f}, 0));
    add_limited_hinge(w, d.torso, d.head, neck, -0.6f, 0.6f, &d.hinges);

    for (int side = 0; side < 2; ++side) {
        // Arms: upper arm hangs from the shoulder, forearm from the elbow (which bends forward).
        const Vec2 shoulder = torso + Vec2{0, 0.45f};
        d.upper_arm[side] = part(make_box({0.1f, 0.28f}, shoulder + Vec2{0, -0.28f}, 0));
        add_limited_hinge(w, d.torso, d.upper_arm[side], shoulder, -2.5f, 2.5f, &d.hinges);
        const Vec2 elbow = shoulder + Vec2{0, -0.56f};
        d.forearm[side] = part(make_box({0.09f, 0.26f}, elbow + Vec2{0, -0.26f}, 0));
        add_limited_hinge(w, d.upper_arm[side], d.forearm[side], elbow, 0.0f, 2.2f, &d.hinges);

        // Legs: thigh hangs from the hip, shin from the knee (which bends backward).
        const Vec2 hip = torso + Vec2{0, -0.5f};
        d.thigh[side] = part(make_box({0.13f, 0.32f}, hip + Vec2{0, -0.32f}, 0));
        add_limited_hinge(w, d.torso, d.thigh[side], hip, -1.2f, 1.2f, &d.hinges);
        const Vec2 knee = hip + Vec2{0, -0.64f};
        d.shin[side] = part(make_box({0.11f, 0.3f}, knee + Vec2{0, -0.3f}, 0));
        add_limited_hinge(w, d.thigh[side], d.shin[side], knee, -2.2f, 0.0f, &d.hinges);
        d.knee_hinge[side] = d.hinges.back();
    }
    return d;
}

Car build_car(World& w, Vec2 pos, int group) {
    Car car;
    Body chassis = make_box({1.3f, 0.28f}, pos, 0, 1.5f);
    set_group(chassis, group);
    chassis.friction = 0.4f;
    car.chassis = w.add(chassis);

    const Real offsets[2] = {-0.95f, 0.95f};  // rear, front
    for (int i = 0; i < 2; ++i) {
        const Vec2 hub = pos + Vec2{offsets[i], -0.62f};

        // The carrier collides with nothing: it only joins the slide, the spring and the wheel's hinge.
        // It is given a wheel-like mass (about 1 kg) on purpose. A spring-damper joint's stiffness is
        // frequency x the effective mass of the two bodies it joins, here chassis and carrier; with a
        // feather-light carrier the spring would be feather-weak and the suspension would bottom out.
        Body carrier = make_box({0.06f, 0.06f}, hub, 0, 80);
        carrier.category = static_cast<std::uint16_t>(1u << group);
        carrier.mask = 0;
        car.carrier[i] = w.add(carrier);

        Body wheel(Shape::make_circle(kCarWheelRadius), hub, 0, BodyType::Dynamic, 1.5f);
        set_group(wheel, group);
        wheel.friction = 1.2f;
        wheel.restitution = 0.05f;
        car.wheel[i] = w.add(wheel);

        // Suspension: slide up and down 25 cm either way, held at rest height by a 2 Hz spring.
        Joint slide = Joint::prismatic(w.bodies, car.chassis, car.carrier[i], hub, {0, 1});
        slide.enable_limit = true;
        slide.lower = -0.25f;
        slide.upper = 0.25f;
        w.add_joint(slide);
        Joint spring = Joint::distance(w.bodies, car.chassis, car.carrier[i], {hub.x, pos.y}, hub);
        spring.frequency_hz = 2.0f;
        spring.damping_ratio = 0.5f;
        w.add_joint(spring);

        // The wheel spins on the carrier; the motor starts out idle.
        Joint axle = Joint::revolute(w.bodies, car.carrier[i], car.wheel[i], hub);
        axle.enable_motor = true;
        axle.motor_speed = 0;
        axle.max_motor = 0;
        car.drive_hinge[i] = w.add_joint(axle);
    }
    return car;
}

void set_car_throttle(World& w, const Car& car, Real wheel_speed, Real max_torque) {
    for (int i = 0; i < 2; ++i) {
        Joint& j = w.joints[static_cast<size_t>(car.drive_hinge[i])];
        j.motor_speed = -wheel_speed;  // a wheel rolling to the right turns clockwise, i.e. negatively
        j.max_motor = max_torque;
    }
}

}  // namespace phys::scenes
