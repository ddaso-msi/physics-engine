#include <phys/joints.hpp>

#include <algorithm>
#include <limits>

namespace phys {

namespace {

constexpr Real kInf = std::numeric_limits<Real>::infinity();

// Largest position error a single correction pass will remove, so a badly violated joint is eased
// back rather than teleported.
constexpr Real kMaxLinearCorrection = static_cast<Real>(0.2);
constexpr Real kMaxAngularCorrection = static_cast<Real>(0.14);  // about 8 degrees

// Angle difference mapped into (-pi, pi], so a joint a hair past +pi does not read as -pi.
Real wrap_angle(Real a) { return a - 2 * kPi * std::round(a / (2 * kPi)); }

Vec2 to_local(const std::vector<Body>& bodies, int index, Vec2 world) {
    return index < 0 ? world : apply_inv(bodies[static_cast<size_t>(index)].transform(), world);
}
Real angle_of(const std::vector<Body>& bodies, int index) {
    return index < 0 ? 0 : bodies[static_cast<size_t>(index)].angle;
}

}  // namespace

// ---- factories -----------------------------------------------------------------------------------

Joint Joint::distance(const std::vector<Body>& bodies, int a, int b, Vec2 world_a, Vec2 world_b) {
    Joint j;
    j.type = JointType::Distance;
    j.a = a;
    j.b = b;
    j.local_anchor_a = to_local(bodies, a, world_a);
    j.local_anchor_b = to_local(bodies, b, world_b);
    j.length = phys::distance(world_a, world_b);
    return j;
}

Joint Joint::revolute(const std::vector<Body>& bodies, int a, int b, Vec2 world_anchor) {
    Joint j;
    j.type = JointType::Revolute;
    j.a = a;
    j.b = b;
    j.local_anchor_a = to_local(bodies, a, world_anchor);
    j.local_anchor_b = to_local(bodies, b, world_anchor);
    j.reference_angle = angle_of(bodies, b) - angle_of(bodies, a);
    return j;
}

Joint Joint::prismatic(const std::vector<Body>& bodies, int a, int b, Vec2 world_anchor, Vec2 world_axis) {
    Joint j;
    j.type = JointType::Prismatic;
    j.a = a;
    j.b = b;
    j.local_anchor_a = to_local(bodies, a, world_anchor);
    j.local_anchor_b = to_local(bodies, b, world_anchor);
    j.reference_angle = angle_of(bodies, b) - angle_of(bodies, a);
    const Vec2 axis = world_axis.normalized();
    j.local_axis_a = a < 0 ? axis : inv_rotate(bodies[static_cast<size_t>(a)].q, axis);
    return j;
}

Joint Joint::mouse(const std::vector<Body>& bodies, int body, Vec2 world_point) {
    Joint j;
    j.type = JointType::Mouse;
    j.a = body;
    j.b = -1;
    j.local_anchor_a = to_local(bodies, body, world_point);
    j.local_anchor_b = world_point;
    j.target = world_point;
    j.frequency_hz = 5;
    j.damping_ratio = static_cast<Real>(0.7);
    j.max_force = 20 * bodies[static_cast<size_t>(body)].mass * static_cast<Real>(9.81);  // 20x its weight
    return j;
}

Vec2 Joint::world_anchor_a(const std::vector<Body>& bodies) const {
    return a < 0 ? local_anchor_a : apply(bodies[static_cast<size_t>(a)].transform(), local_anchor_a);
}
Vec2 Joint::world_anchor_b(const std::vector<Body>& bodies) const {
    return b < 0 ? local_anchor_b : apply(bodies[static_cast<size_t>(b)].transform(), local_anchor_b);
}

// ---- rows ----------------------------------------------------------------------------------------

namespace {

// The velocity change a row's impulse causes: M^-1 J^T lambda, applied to both bodies.
void apply_impulse(JointRow& r, Real lambda) {
    r.a->vel += r.lin_a * (r.a->inv_mass * lambda);
    r.a->w += r.a->inv_inertia * r.ang_a * lambda;
    r.b->vel += r.lin_b * (r.b->inv_mass * lambda);
    r.b->w += r.b->inv_inertia * r.ang_b * lambda;
}

// The same push applied to POSITIONS instead: M^-1 J^T lambda read as a displacement. No velocity
// changes, so correcting an error this way cannot inject energy.
void apply_displacement(JointRow& r, Real lambda) {
    r.a->pos += r.lin_a * (r.a->inv_mass * lambda);
    r.a->set_angle(r.a->angle + r.a->inv_inertia * r.ang_a * lambda);
    r.b->pos += r.lin_b * (r.b->inv_mass * lambda);
    r.b->set_angle(r.b->angle + r.b->inv_inertia * r.ang_b * lambda);
}

}  // namespace

void JointSolver::prepare(std::vector<Body>& bodies, std::span<Joint> joints, Real dt, const SolverSettings& s,
                          bool for_velocity) {
    rows_.clear();
    world_ = Body();

    for (Joint& j : joints) {
        Body* A = j.a < 0 ? &world_ : &bodies[static_cast<size_t>(j.a)];
        Body* B = j.b < 0 ? &world_ : &bodies[static_cast<size_t>(j.b)];

        // Velocity pass: move last step's impulses aside; slots not claimed by a row this step end up
        // zero. (The position pass must leave the stored impulses alone.)
        Real old[Joint::kMaxRows] = {};
        if (for_velocity) {
            std::copy(std::begin(j.impulses), std::end(j.impulses), old);
            std::fill(std::begin(j.impulses), std::end(j.impulses), Real(0));
        }

        // Adds one row. `bias` is the extra velocity the row should leave (a motor's target speed, or a
        // limit's allowance). `error` is the position error this row measures, for the position pass.
        auto add = [&](int slot, Vec2 lin_a, Real ang_a, Vec2 lin_b, Real ang_b, Real bias, Real lower = -kInf,
                       Real upper = kInf, Real error = 0) -> JointRow& {
            JointRow r;
            r.a = A;
            r.b = B;
            r.lin_a = lin_a;
            r.ang_a = ang_a;
            r.lin_b = lin_b;
            r.ang_b = ang_b;
            r.bias = bias;
            r.error = error;
            r.lower = lower;
            r.upper = upper;
            const Real k = A->inv_mass * dot(lin_a, lin_a) + A->inv_inertia * ang_a * ang_a +
                           B->inv_mass * dot(lin_b, lin_b) + B->inv_inertia * ang_b * ang_b;
            r.mass = k > 0 ? 1 / k : 0;
            if (for_velocity) {
                r.store = &j.impulses[slot];
                r.impulse = s.warm_starting ? std::clamp(old[slot], lower, upper) : Real(0);
            }
            rows_.push_back(r);
            return rows_.back();
        };

        // A spring-damper row. Instead of forcing C to zero this step, treat the row as a spring of
        // the given frequency and damping ratio, using the implicit-Euler form so it stays stable
        // however stiff it is:  the restoring term is  C * dt*k / (d + dt*k)  and the effective mass
        // is softened by gamma = 1 / (dt (d + dt k)), with k = m w^2 and d = 2 m zeta w for the row's
        // own effective mass m. A spring has no "position error" to remove: it is meant to stretch.
        auto add_soft = [&](int slot, Vec2 lin_a, Real ang_a, Vec2 lin_b, Real ang_b, Real C, Real hz, Real zeta,
                            Real lower = -kInf, Real upper = kInf) {
            JointRow& r = add(slot, lin_a, ang_a, lin_b, ang_b, 0, lower, upper);
            if (r.mass <= 0) return;
            const Real m = r.mass, omega = 2 * kPi * hz;
            const Real d = 2 * m * zeta * omega, k = m * omega * omega;
            const Real gamma = 1 / (dt * (d + dt * k));
            r.softness = gamma;
            r.bias = C * dt * k * gamma;
            r.mass = 1 / (1 / m + gamma);
        };

        // Limit rows are one-sided. While the limit is not reached (C > 0) the bodies may close in on it,
        // but no faster than would reach it this step. Once violated (C < 0) the velocity pass only
        // forbids moving further in; the position pass is what brings the bodies back inside.
        auto limit_bias = [&](Real C) { return C > 0 ? C / dt : Real(0); };

        const Vec2 ra = rotate(A->q, j.local_anchor_a), rb = rotate(B->q, j.local_anchor_b);
        const Vec2 pa = A->pos + ra, pb = B->pos + rb;

        switch (j.type) {
            case JointType::Distance: {
                const Vec2 d = pb - pa;
                const Real len = d.length();
                const Vec2 u = len > static_cast<Real>(1e-6) ? d * (1 / len) : Vec2{1, 0};
                const Real C = len - j.length;
                const Vec2 lin_a = -u, lin_b = u;
                const Real ang_a = -cross(ra, u), ang_b = cross(rb, u);
                if (j.frequency_hz > 0) add_soft(0, lin_a, ang_a, lin_b, ang_b, C, j.frequency_hz, j.damping_ratio);
                else add(0, lin_a, ang_a, lin_b, ang_b, 0, -kInf, kInf, C);
                break;
            }

            case JointType::Revolute: {
                // Keep the two anchor points together: one row per world axis.
                const Vec2 C = pb - pa;
                const Vec2 axes[2] = {{1, 0}, {0, 1}};
                for (int i = 0; i < 2; ++i) {
                    const Vec2 e = axes[i];
                    add(i, -e, -cross(ra, e), e, cross(rb, e), 0, -kInf, kInf, dot(C, e));
                }
                if (j.enable_motor)
                    add(2, {}, -1, {}, 1, -j.motor_speed, -j.max_motor * dt, j.max_motor * dt);
                if (j.enable_limit) {
                    const Real angle = wrap_angle(B->angle - A->angle - j.reference_angle);
                    const Real c_low = angle - j.lower, c_up = j.upper - angle;
                    add(3, {}, -1, {}, 1, limit_bias(c_low), 0, kInf, std::min(c_low, Real(0)));
                    add(4, {}, 1, {}, -1, limit_bias(c_up), 0, kInf, std::min(c_up, Real(0)));
                }
                break;
            }

            case JointType::Prismatic: {
                const Vec2 d = pb - pa;
                const Vec2 axis = rotate(A->q, j.local_axis_a);
                const Vec2 perp = phys::perp(axis);

                // Stay on the rail: no motion of the anchors perpendicular to the axis. The axis turns
                // with body a, so the row also has a term for that turning (the s1 below).
                {
                    const Real s1 = cross(d + ra, perp), s2 = cross(rb, perp);
                    add(0, -perp, -s1, perp, s2, 0, -kInf, kInf, dot(perp, d));
                }
                // No relative rotation.
                add(1, {}, -1, {}, 1, 0, -kInf, kInf, wrap_angle(B->angle - A->angle - j.reference_angle));

                // Along the rail.
                const Real a1 = cross(d + ra, axis), a2 = cross(rb, axis);
                const Real t = dot(axis, d);
                if (j.enable_motor)
                    add(2, -axis, -a1, axis, a2, -j.motor_speed, -j.max_motor * dt, j.max_motor * dt);
                if (j.enable_limit) {
                    const Real c_low = t - j.lower, c_up = j.upper - t;
                    add(3, -axis, -a1, axis, a2, limit_bias(c_low), 0, kInf, std::min(c_low, Real(0)));
                    add(4, axis, a1, -axis, -a2, limit_bias(c_up), 0, kInf, std::min(c_up, Real(0)));
                }
                break;
            }

            case JointType::Mouse: {
                // The grabbed point of body a is pulled toward `target`; b is the world. Impulses are
                // capped so a heavy body is dragged, not flung. (A spring: no position error to remove.)
                const Vec2 C = pa - j.target;
                const Vec2 axes[2] = {{1, 0}, {0, 1}};
                const Real cap = j.max_force * dt;
                for (int i = 0; i < 2; ++i) {
                    const Vec2 e = axes[i];
                    if (j.frequency_hz > 0)
                        add_soft(i, e, cross(ra, e), {}, 0, dot(C, e), j.frequency_hz, j.damping_ratio, -cap, cap);
                    else
                        add(i, e, cross(ra, e), {}, 0, 0, -cap, cap, dot(C, e));
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
        const Real jv = dot(r.lin_a, r.a->vel) + r.ang_a * r.a->w + dot(r.lin_b, r.b->vel) + r.ang_b * r.b->w;
        Real lambda = -r.mass * (jv + r.bias + r.softness * r.impulse);
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
            // Rows for this one joint, built from the poses as they are right now.
            solver.prepare(bodies, std::span<Joint>(&joint, 1), 1, settings, /*for_velocity=*/false);
            for (JointRow& r : solver.rows_mutable()) {
                if (r.error == 0 || r.mass <= 0) continue;
                const bool angular_only = r.lin_a.length_sq() == 0 && r.lin_b.length_sq() == 0;
                const Real limit = angular_only ? kMaxAngularCorrection : kMaxLinearCorrection;
                const Real C = std::clamp(r.error, -limit, limit);
                apply_displacement(r, -C * r.mass);
            }
        }
    }
}

}  // namespace phys
