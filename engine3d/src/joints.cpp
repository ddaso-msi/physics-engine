#include <phys3d/joints.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace phys3d {

namespace {

constexpr Real kInf = std::numeric_limits<Real>::infinity();

// Largest position error a single correction will remove, so a badly violated joint is eased back rather
// than teleported.
constexpr Real kMaxLinearCorrection = static_cast<Real>(0.2);
constexpr Real kMaxAngularCorrection = static_cast<Real>(0.14);  // about 8 degrees

Vec3 to_local(const std::vector<Body>& bodies, int index, Vec3 world) {
    return index < 0 ? world : apply_inv(bodies[static_cast<size_t>(index)].transform(), world);
}
Vec3 dir_to_local(const std::vector<Body>& bodies, int index, Vec3 world) {
    return index < 0 ? world : inv_rotate(bodies[static_cast<size_t>(index)].q, world);
}
Vec3 dir_to_world(const std::vector<Body>& bodies, int index, Vec3 local) {
    return index < 0 ? local : rotate(bodies[static_cast<size_t>(index)].q, local);
}

// How far `ref_b` has turned from `ref_a` about unit `axis` (both perpendicular to it), in (-pi, pi].
Real angle_about(Vec3 axis, Vec3 ref_a, Vec3 ref_b) { return std::atan2(dot(cross(ref_a, ref_b), axis), dot(ref_a, ref_b)); }

}  // namespace

// ---- factories -----------------------------------------------------------------------------------

Joint Joint::distance(const std::vector<Body>& bodies, int a, int b, Vec3 world_a, Vec3 world_b) {
    Joint j;
    j.type = JointType::Distance;
    j.a = a;
    j.b = b;
    j.local_anchor_a = to_local(bodies, a, world_a);
    j.local_anchor_b = to_local(bodies, b, world_b);
    j.length = phys3d::distance(world_a, world_b);
    return j;
}

Joint Joint::ball(const std::vector<Body>& bodies, int a, int b, Vec3 world_anchor) {
    Joint j;
    j.type = JointType::Ball;
    j.a = a;
    j.b = b;
    j.local_anchor_a = to_local(bodies, a, world_anchor);
    j.local_anchor_b = to_local(bodies, b, world_anchor);
    return j;
}

Joint Joint::hinge(const std::vector<Body>& bodies, int a, int b, Vec3 world_anchor, Vec3 world_axis) {
    Joint j = ball(bodies, a, b, world_anchor);
    j.type = JointType::Hinge;
    const Vec3 axis = world_axis.normalized();
    // Any direction perpendicular to the axis will do as the zero mark: cross the axis with whichever
    // world axis it is least aligned with.
    const Vec3 ref = cross(axis, std::fabs(axis.x) < static_cast<Real>(0.57735) ? Vec3{1, 0, 0} : Vec3{0, 1, 0}).normalized();
    j.local_axis_a = dir_to_local(bodies, a, axis);
    j.local_axis_b = dir_to_local(bodies, b, axis);
    j.local_ref_a = dir_to_local(bodies, a, ref);
    j.local_ref_b = dir_to_local(bodies, b, ref);
    return j;
}

Vec3 Joint::world_anchor_a(const std::vector<Body>& bodies) const {
    return a < 0 ? local_anchor_a : apply(bodies[static_cast<size_t>(a)].transform(), local_anchor_a);
}
Vec3 Joint::world_anchor_b(const std::vector<Body>& bodies) const {
    return b < 0 ? local_anchor_b : apply(bodies[static_cast<size_t>(b)].transform(), local_anchor_b);
}
Vec3 Joint::world_axis(const std::vector<Body>& bodies) const { return dir_to_world(bodies, a, local_axis_a); }
Real Joint::hinge_angle(const std::vector<Body>& bodies) const {
    return angle_about(world_axis(bodies), dir_to_world(bodies, a, local_ref_a), dir_to_world(bodies, b, local_ref_b));
}

// ---- rows ----------------------------------------------------------------------------------------

namespace {

// The velocity change a row's impulse causes: M^-1 J^T lambda, applied to both bodies.
void apply_impulse(JointRow& r, Real lambda) {
    r.a->vel += r.lin_a * (r.a->inv_mass * lambda);
    r.a->w += r.inv_i_ang_a * lambda;
    r.b->vel += r.lin_b * (r.b->inv_mass * lambda);
    r.b->w += r.inv_i_ang_b * lambda;
}

// The same push applied to POSES instead: M^-1 J^T lambda read as a displacement and a small rotation.
void apply_displacement(JointRow& r, Real lambda) {
    r.a->pos += r.lin_a * (r.a->inv_mass * lambda);
    r.a->q = integrate_orientation(r.a->q, r.inv_i_ang_a * lambda, 1);
    r.b->pos += r.lin_b * (r.b->inv_mass * lambda);
    r.b->q = integrate_orientation(r.b->q, r.inv_i_ang_b * lambda, 1);
}

}  // namespace

void JointSolver::prepare(std::vector<Body>& bodies, std::span<Joint> joints, Real dt, const SolverSettings& s,
                          bool for_velocity) {
    rows_.clear();
    world_ = Body::fixed({});

    for (Joint& j : joints) {
        // A joint whose bodies are all asleep (or fixed) has nothing to do; its stored impulses wait.
        const bool live = (j.a >= 0 && bodies[static_cast<size_t>(j.a)].is_active()) ||
                          (j.b >= 0 && bodies[static_cast<size_t>(j.b)].is_active());
        if (!live) continue;
        Body* A = j.a < 0 ? &world_ : &bodies[static_cast<size_t>(j.a)];
        Body* B = j.b < 0 ? &world_ : &bodies[static_cast<size_t>(j.b)];
        const Mat3 inv_ia = A->inv_inertia_world(), inv_ib = B->inv_inertia_world();

        // Velocity pass: move last step's impulses aside; slots not claimed by a row this step end up zero.
        Real old[Joint::kMaxRows] = {};
        if (for_velocity) {
            std::copy(std::begin(j.impulses), std::end(j.impulses), old);
            std::fill(std::begin(j.impulses), std::end(j.impulses), Real(0));
        }

        // Adds one row. `bias` is extra velocity the row should leave (a motor's speed, a limit's
        // allowance); `error` is the position error it measures, for the position pass.
        auto add = [&](int slot, Vec3 lin_a, Vec3 ang_a, Vec3 lin_b, Vec3 ang_b, Real bias, Real lower = -kInf,
                       Real upper = kInf, Real error = 0) -> JointRow& {
            JointRow r;
            r.a = A;
            r.b = B;
            r.lin_a = lin_a;
            r.ang_a = ang_a;
            r.lin_b = lin_b;
            r.ang_b = ang_b;
            r.inv_i_ang_a = inv_ia * ang_a;
            r.inv_i_ang_b = inv_ib * ang_b;
            r.bias = bias;
            r.error = error;
            r.lower = lower;
            r.upper = upper;
            const Real k = A->inv_mass * dot(lin_a, lin_a) + dot(ang_a, r.inv_i_ang_a) +
                           B->inv_mass * dot(lin_b, lin_b) + dot(ang_b, r.inv_i_ang_b);
            r.mass = k > 0 ? 1 / k : 0;
            if (for_velocity) {
                r.store = &j.impulses[slot];
                r.impulse = s.warm_starting ? std::clamp(old[slot], lower, upper) : Real(0);
            }
            rows_.push_back(r);
            return rows_.back();
        };

        const Vec3 ra = rotate(A->q, j.local_anchor_a), rb = rotate(B->q, j.local_anchor_b);
        const Vec3 pa = A->pos + ra, pb = B->pos + rb;

        // Keeps the two anchor points together along world direction e. The anchor of b moves at
        // v_b + w_b x rb, and (w x r) . e = w . (r x e), which gives the angular parts.
        auto add_point_row = [&](int slot, Vec3 e, Real error) {
            add(slot, -e, -cross(ra, e), e, cross(rb, e), 0, -kInf, kInf, error);
        };

        switch (j.type) {
            case JointType::Distance: {
                const Vec3 d = pb - pa;
                const Real len = d.length();
                const Vec3 u = len > static_cast<Real>(1e-6) ? d * (1 / len) : Vec3{1, 0, 0};
                const Real C = len - j.length;
                if (j.frequency_hz > 0) {
                    // A spring-damper in implicit-Euler form, as in 2D: stable however stiff. A spring has
                    // no position error to remove; it is meant to stretch.
                    JointRow& r = add(0, -u, -cross(ra, u), u, cross(rb, u), 0);
                    if (r.mass > 0) {
                        const Real m = r.mass, omega = 2 * kPi * j.frequency_hz;
                        const Real damping = 2 * m * j.damping_ratio * omega, stiffness = m * omega * omega;
                        const Real gamma = 1 / (dt * (damping + dt * stiffness));
                        r.softness = gamma;
                        r.bias = C * dt * stiffness * gamma;
                        r.mass = 1 / (1 / m + gamma);
                    }
                } else {
                    add_point_row(0, u, C);
                }
                break;
            }

            case JointType::Ball:
            case JointType::Hinge: {
                // The point: one row per world axis.
                const Vec3 C = pb - pa;
                const Vec3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
                for (int i = 0; i < 3; ++i) add_point_row(i, axes[i], dot(C, axes[i]));
                if (j.type == JointType::Ball) break;

                // The axis. With a1 and b1 the axis as each body carries it, a1 x b1 is (for a small
                // misalignment) the rotation that has pulled them apart, and it is perpendicular to a1.
                // So: no relative turning about the two directions perpendicular to a1. Those two
                // directions ride with body a, so the rows keep their meaning from step to step.
                const Vec3 a1 = rotate(A->q, j.local_axis_a), b1 = rotate(B->q, j.local_axis_b);
                const Vec3 t1 = rotate(A->q, j.local_ref_a), t2 = cross(a1, t1);
                const Vec3 apart = cross(a1, b1);
                add(3, {}, -t1, {}, t1, 0, -kInf, kInf, dot(apart, t1));
                add(4, {}, -t2, {}, t2, 0, -kInf, kInf, dot(apart, t2));

                // About the axis: a motor asks for a relative speed and may use up to max_motor of torque.
                if (j.enable_motor) add(5, {}, -a1, {}, a1, -j.motor_speed, -j.max_motor * dt, j.max_motor * dt);
                if (j.enable_limit) {
                    // One-sided rows. Short of the limit (C > 0) the bodies may close in on it, but no
                    // faster than would reach it this step. Past it (C < 0) the velocity pass only forbids
                    // going further; the position pass brings them back.
                    const Real angle = angle_about(a1, t1, rotate(B->q, j.local_ref_b));
                    const Real c_low = angle - j.lower, c_up = j.upper - angle;
                    add(6, {}, -a1, {}, a1, c_low > 0 ? c_low / dt : Real(0), 0, kInf, std::min(c_low, Real(0)));
                    add(7, {}, a1, {}, -a1, c_up > 0 ? c_up / dt : Real(0), 0, kInf, std::min(c_up, Real(0)));
                }
                break;
            }
        }
    }
}

void JointSolver::warm_start() {
    for (JointRow& r : rows_) apply_impulse(r, r.impulse);
}

void JointSolver::solve_velocity() {
    for (JointRow& r : rows_) {
        Real lambda = -r.mass * (r.velocity() + r.bias + r.softness * r.impulse);
        const Real total = std::clamp(r.impulse + lambda, r.lower, r.upper);  // clamp the ACCUMULATED impulse
        lambda = total - r.impulse;
        r.impulse = total;
        apply_impulse(r, lambda);
    }
}

void JointSolver::store() {
    for (const JointRow& r : rows_) *r.store = r.impulse;
}

void solve_joint_positions(std::vector<Body>& bodies, std::vector<Joint>& joints, const SolverSettings& settings) {
    JointSolver solver;
    for (int pass = 0; pass < settings.joint_position_iterations; ++pass) {
        for (Joint& joint : joints) {
            // One ROW at a time, each measured from the poses as they are right now. A joint's rows share
            // its bodies, so correcting one changes the error of the others: closing a hinge's point turns
            // the body, which changes its angle, and turning it back to a limit moves the point. Corrected
            // all together from one measurement they overshoot against each other and never settle.
            for (size_t k = 0;; ++k) {
                solver.prepare(bodies, std::span<Joint>(&joint, 1), 1, settings, /*for_velocity=*/false);
                if (k >= solver.rows().size()) break;
                JointRow& r = solver.rows_mutable()[k];
                if (r.error == 0 || r.mass <= 0) continue;
                const bool angular_only = r.lin_a.length_sq() == 0 && r.lin_b.length_sq() == 0;
                const Real limit = angular_only ? kMaxAngularCorrection : kMaxLinearCorrection;
                const Real C = std::clamp(r.error, -limit, limit);
                apply_displacement(r, -C * r.mass);
            }
        }
    }
}

}  // namespace phys3d
