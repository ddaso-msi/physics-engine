#include <phys3d/joints.hpp>

#include "joint_trees.hpp"
#include "solving_order.hpp"

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
// Errors below which the joint trees' position pass is not repeated: a millimetre, or a milliradian.
constexpr Real kTreesSettled = static_cast<Real>(1e-3);

Vec3 to_local(const std::vector<Body>& bodies, int index, Vec3 world) {
    return index < 0 ? world : apply_inv(bodies[static_cast<size_t>(index)].transform(), world);
}
Vec3 dir_to_local(const std::vector<Body>& bodies, int index, Vec3 world) {
    return index < 0 ? world : inv_rotate(bodies[static_cast<size_t>(index)].q, world);
}
Vec3 dir_to_world(const std::vector<Body>& bodies, int index, Vec3 local) {
    return index < 0 ? local : rotate(bodies[static_cast<size_t>(index)].q, local);
}

Quat orientation_of(const std::vector<Body>& bodies, int index) { return index < 0 ? Quat{} : bodies[static_cast<size_t>(index)].q; }

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

Joint Joint::weld(const std::vector<Body>& bodies, int a, int b, Vec3 world_anchor) {
    Joint j = ball(bodies, a, b, world_anchor);
    j.type = JointType::Weld;
    j.relative_orientation = conjugate(orientation_of(bodies, a)) * orientation_of(bodies, b);
    return j;
}

Joint Joint::slider(const std::vector<Body>& bodies, int a, int b, Vec3 world_anchor, Vec3 world_axis) {
    Joint j = hinge(bodies, a, b, world_anchor, world_axis);  // the same anchors, axis and reference directions
    j.type = JointType::Slider;
    j.relative_orientation = conjugate(orientation_of(bodies, a)) * orientation_of(bodies, b);
    return j;
}

void Joint::limit_cone(const std::vector<Body>& bodies, Vec3 world_axis, Real max_angle) {
    const Vec3 axis = world_axis.normalized();
    local_axis_a = dir_to_local(bodies, a, axis);
    local_axis_b = dir_to_local(bodies, b, axis);
    enable_cone = true;
    cone_angle = max_angle;
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

Real Joint::slide(const std::vector<Body>& bodies) const { return dot(world_anchor_b(bodies) - world_anchor_a(bodies), world_axis(bodies)); }
Real Joint::cone_tilt(const std::vector<Body>& bodies) const {
    const Vec3 a1 = world_axis(bodies), b1 = dir_to_world(bodies, b, local_axis_b);
    return std::atan2(cross(a1, b1).length(), dot(a1, b1));
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

// A joint takes part in solving unless all its bodies are asleep or fixed.
bool is_live(const std::vector<Body>& bodies, const Joint& j) {
    return (j.a >= 0 && bodies[static_cast<size_t>(j.a)].is_active()) || (j.b >= 0 && bodies[static_cast<size_t>(j.b)].is_active());
}

// The position pass for a joint's anchor point (ball, hinge and weld): the three rows that hold the two anchors
// together, corrected as one. The gap C = pb - pa is measured once, and the push P that closes it is found
// from the 3 x 3 matrix K that says how the gap answers a push:
//     a push P on b at its anchor (and -P on a) moves the gap by  K P,
//     K = (1/ma + 1/mb) E - [ra] IA^-1 [ra] - [rb] IB^-1 [rb]       ([r] the matrix of r x)
// so P = -K^-1 C. The three rows one at a time, each measured afresh, come to the same answer only in the limit:
// each sees the others' rows only through what they have already moved. (Step 21.)
void correct_anchor(std::vector<Body>& bodies, const Joint& j) {
    Body* A = j.a < 0 ? nullptr : &bodies[static_cast<size_t>(j.a)];
    Body* B = j.b < 0 ? nullptr : &bodies[static_cast<size_t>(j.b)];
    const Vec3 ra = A ? rotate(A->q, j.local_anchor_a) : Vec3{}, rb = B ? rotate(B->q, j.local_anchor_b) : Vec3{};
    const Vec3 pa = A ? A->pos + ra : j.local_anchor_a, pb = B ? B->pos + rb : j.local_anchor_b;
    Vec3 C = pb - pa;
    const Real gap = C.length();
    if (gap == 0) return;
    if (gap > kMaxLinearCorrection) C = C * (kMaxLinearCorrection / gap);
    const bool a_moves = A && A->type == BodyType::Dynamic, b_moves = B && B->type == BodyType::Dynamic;
    const Mat3 inv_ia = a_moves ? A->inv_inertia_world() : Mat3::diagonal(0, 0, 0);
    const Mat3 inv_ib = b_moves ? B->inv_inertia_world() : Mat3::diagonal(0, 0, 0);
    const Real ma = a_moves ? A->inv_mass : 0, mb = b_moves ? B->inv_mass : 0;
    const Mat3 K = Mat3::identity() * (ma + mb) - Mat3::skew(ra) * inv_ia * Mat3::skew(ra) - Mat3::skew(rb) * inv_ib * Mat3::skew(rb);
    const Vec3 P = K.inverse() * -C;
    if (a_moves) {
        A->pos -= P * ma;
        A->q = integrate_orientation(A->q, -(inv_ia * cross(ra, P)), 1);
    }
    if (b_moves) {
        B->pos += P * mb;
        B->q = integrate_orientation(B->q, inv_ib * cross(rb, P), 1);
    }
}

// The same push applied to POSES instead: M^-1 J^T lambda read as a displacement and a small rotation. (Not
// to a body that cannot move: for one, both are zero, and turning it by nothing would only renormalise it.)
void apply_displacement(JointRow& r, Real lambda) {
    if (r.a->type == BodyType::Dynamic) {
        r.a->pos += r.lin_a * (r.a->inv_mass * lambda);
        r.a->q = integrate_orientation(r.a->q, r.inv_i_ang_a * lambda, 1);
    }
    if (r.b->type == BodyType::Dynamic) {
        r.b->pos += r.lin_b * (r.b->inv_mass * lambda);
        r.b->q = integrate_orientation(r.b->q, r.inv_i_ang_b * lambda, 1);
    }
}

}  // namespace

JointSolver::JointSolver() : trees_(std::make_unique<JointTrees>()) {}
JointSolver::~JointSolver() = default;
JointSolver::JointSolver(JointSolver&&) noexcept = default;
JointSolver& JointSolver::operator=(JointSolver&&) noexcept = default;

void JointSolver::prepare(std::vector<Body>& bodies, std::span<Joint> joints, Real dt, const SolverSettings& s,
                          bool for_velocity, int only_row) {
    rows_.clear();
    counted_ = 0;
    first_row_.clear();
    joint_of_.clear();

    // A joint whose bodies are all asleep (or fixed) has nothing to do; its stored impulses wait. The rest are
    // taken colour by colour (see solving_order.hpp), unless there is only one: the position pass prepares
    // one joint at a time, many times a step, and has no order to choose.
    auto live = [&](const Joint& j) { return is_live(bodies, j); };
    const bool coloured = joints.size() > 1;
    std::vector<size_t> order;
    if (coloured) order = solving_order(bodies, joints.size(), [&](size_t i, int& a, int& b) {
        a = joints[i].a, b = joints[i].b;
        return live(joints[i]);
    });
    for (size_t k = 0; k < (coloured ? order.size() : joints.size()); ++k) {
        Joint& j = joints[coloured ? order[k] : k];
        if (!coloured && !live(j)) continue;
        Body* A = j.a < 0 ? &world_ : &bodies[static_cast<size_t>(j.a)];
        Body* B = j.b < 0 ? &world_ : &bodies[static_cast<size_t>(j.b)];
        const Mat3 inv_ia = A->inv_inertia_world(), inv_ib = B->inv_inertia_world();
        first_row_.push_back(static_cast<int>(rows_.size()));
        joint_of_.push_back(static_cast<int>(&j - joints.data()));

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
            if (only_row >= 0 && counted_++ != only_row) return skipped_;
            if (only_row < 0) ++counted_;
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
                if (j.type == JointType::Ball) {
                    if (!j.enable_cone) break;
                    // The cone: the angle between the two copies of the axis may not exceed cone_angle. That
                    // angle grows when b turns, relative to a, about the direction a1 x b1; so this is a
                    // one-sided row about that direction, like a hinge's limit. (With the two copies in
                    // line the direction is undefined, and the limit is nowhere near.)
                    const Vec3 a1 = rotate(A->q, j.local_axis_a), b1 = rotate(B->q, j.local_axis_b);
                    const Vec3 across = cross(a1, b1);
                    const Real sine = across.length();
                    if (sine < static_cast<Real>(1e-6)) break;
                    const Vec3 u = across / sine;
                    const Real room = j.cone_angle - std::atan2(sine, dot(a1, b1));
                    add(3, {}, u, {}, -u, room > 0 ? room / dt : Real(0), 0, kInf, std::min(room, Real(0)));
                    break;
                }

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

            case JointType::Weld:
            case JointType::Slider: {
                // No relative turning. Seen from a, b's orientation should be what it was when the joint
                // was made; whatever rotation takes that to what it is now is the error, a small turn
                // about some axis, and each world axis gets a row for its share of it.
                const Quat seen_from_a = conjugate(A->q) * B->q;
                Quat off = seen_from_a * conjugate(j.relative_orientation);
                if (off.w < 0) off = off * Real(-1);  // the same rotation, written the short way round
                const Vec3 v = off.vector_part();
                const Real sine = v.length();
                const Vec3 turned = sine > static_cast<Real>(1e-9) ? rotate(A->q, v * (2 * std::atan2(sine, off.w) / sine)) : Vec3{};
                const Vec3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
                const int first_angular = j.type == JointType::Weld ? 3 : 2;
                for (int i = 0; i < 3; ++i) add(first_angular + i, {}, -axes[i], {}, axes[i], 0, -kInf, kInf, dot(turned, axes[i]));

                if (j.type == JointType::Weld) {
                    const Vec3 C = pb - pa;
                    for (int i = 0; i < 3; ++i) add_point_row(i, axes[i], dot(C, axes[i]));
                    break;
                }

                // The slider. Its axis and the two directions across it ride with body a. A row along
                // direction e of that frame measures C = e . (pb - pa), and since e itself turns with a,
                //     dC/dt = e . (v_b + w_b x rb - v_a - w_a x ra) + (w_a x e) . d
                // which puts ra + d, not ra, in a's angular part.
                const Vec3 d = pb - pa;
                const Vec3 a1 = rotate(A->q, j.local_axis_a), t1 = rotate(A->q, j.local_ref_a), t2 = cross(a1, t1);
                auto add_rail_row = [&](int slot, Vec3 e, Real bias, Real lower, Real upper, Real error) {
                    add(slot, -e, -cross(ra + d, e), e, cross(rb, e), bias, lower, upper, error);
                };
                // Across the axis: stay on the rail.
                add_rail_row(0, t1, 0, -kInf, kInf, dot(d, t1));
                add_rail_row(1, t2, 0, -kInf, kInf, dot(d, t2));
                // Along it: a motor (a speed, with a force limit) and the ends of the rail.
                if (j.enable_motor) add_rail_row(5, a1, -j.motor_speed, -j.max_motor * dt, j.max_motor * dt, 0);
                if (j.enable_limit) {
                    const Real slid = dot(d, a1), c_low = slid - j.lower, c_up = j.upper - slid;
                    add_rail_row(6, a1, c_low > 0 ? c_low / dt : Real(0), 0, kInf, std::min(c_low, Real(0)));
                    add_rail_row(7, -a1, c_up > 0 ? c_up / dt : Real(0), 0, kInf, std::min(c_up, Real(0)));
                }
                break;
            }
        }
    }
    if (for_velocity) {
        first_row_.push_back(static_cast<int>(rows_.size()));
        trees_->build(bodies, rows_, first_row_);
    }
}

void JointSolver::warm_start() {
    for (JointRow& r : rows_) apply_impulse(r, r.impulse);
}

void JointSolver::solve_velocity(bool exactly) {
    // Every row, one at a time. A joint's rows come one after another and all act on the same two bodies, so
    // their velocities are read once, the rows solved on copies held in registers, and written back once, and
    // only to a body that can move (as the contacts do; see solver.cpp).
    for (size_t first = 0; first < rows_.size();) {
        Body* a = rows_[first].a;
        Body* b = rows_[first].b;
        Vec3 va = a->vel, wa = a->w, vb = b->vel, wb = b->w;
        size_t i = first;
        for (; i < rows_.size() && rows_[i].a == a && rows_[i].b == b; ++i) {
            JointRow& r = rows_[i];
            const Real jv = dot(r.lin_a, va) + dot(r.ang_a, wa) + dot(r.lin_b, vb) + dot(r.ang_b, wb);
            Real lambda = -r.mass * (jv + r.bias + r.softness * r.impulse);
            const Real total = std::clamp(r.impulse + lambda, r.lower, r.upper);  // clamp the ACCUMULATED impulse
            lambda = total - r.impulse;
            r.impulse = total;
            va += r.lin_a * (a->inv_mass * lambda);
            wa += r.inv_i_ang_a * lambda;
            vb += r.lin_b * (b->inv_mass * lambda);
            wb += r.inv_i_ang_b * lambda;
        }
        if (a->type == BodyType::Dynamic) a->vel = va, a->w = wa;
        if (b->type == BodyType::Dynamic) b->vel = vb, b->w = wb;
        first = i;
    }
    // Then the trees, exactly. (Every sweep cost 0.2 ms more at a thousand joints and held them no better: the
    // mean gap in the chains scene was the same, the worst moved both ways by a millimetre or so.)
    if (exactly) trees_->solve(rows_);
}

void JointSolver::correct_trees(std::vector<Body>& bodies, Real max_linear, Real max_angular, std::vector<std::uint8_t>& taken,
                                std::vector<Real>& error) {
    std::vector<int> first = first_row_;
    first.push_back(static_cast<int>(rows_.size()));
    std::vector<std::uint8_t> by_position;
    trees_->build(bodies, rows_, first, &by_position);
    std::vector<Real> largest(by_position.size(), 0);
    trees_->correct(rows_, max_linear, max_angular, largest);
    for (size_t k = 0; k < by_position.size(); ++k)
        if (by_position[k]) {
            taken[static_cast<size_t>(joint_of_[k])] = 1;
            error[static_cast<size_t>(joint_of_[k])] = largest[k];
        }
}

void JointSolver::store() {
    for (const JointRow& r : rows_) *r.store = r.impulse;
}

void solve_joint_positions(std::vector<Body>& bodies, std::vector<Joint>& joints, const SolverSettings& settings) {
    JointSolver solver, all;
    std::vector<std::uint8_t> taken, sub_taken;
    std::vector<Real> error, sub_error;
    std::vector<size_t> again, still;  // the joints of trees that need another pass
    std::vector<Joint> sub;
    for (int pass = 0; pass < settings.joint_position_iterations; ++pass) {
        // First the joints that form trees, all of each tree at once (step 26): their rows that hold both ways
        // measured once and corrected together, so that the correction of one joint is made knowing what it does
        // to every other. (Done joint by joint, step 22 found the opposite: worse, once chains collide.) A
        // correction is exact to first order, so what it leaves is of the order of the square of what it found:
        // after errors under a millimetre (or a milliradian) another pass would find almost nothing, and a tree
        // that had none bigger gets none. (Measured in the chains scene: one pass always left gaps of 6 * 10^-7 m on
        // average, but a few of a millimetre or two where a joint had been far off; four passes for every tree
        // cost twice as much, and so did a second for all of them whenever any one needed it.)
        if (pass == 0) {
            all.prepare(bodies, joints, 1, settings, /*for_velocity=*/false);
            taken.assign(joints.size(), 0);
            error.assign(joints.size(), 0);
            all.correct_trees(bodies, kMaxLinearCorrection, kMaxAngularCorrection, taken, error);
            again.clear();
            for (size_t j = 0; j < joints.size(); ++j)
                if (taken[j] && error[j] > kTreesSettled) again.push_back(j);
        } else if (!again.empty()) {
            sub.clear();
            for (const size_t j : again) sub.push_back(joints[j]);  // (copies: the pass only moves bodies)
            all.prepare(bodies, sub, 1, settings, /*for_velocity=*/false);
            sub_taken.assign(sub.size(), 0);
            sub_error.assign(sub.size(), 0);
            all.correct_trees(bodies, kMaxLinearCorrection, kMaxAngularCorrection, sub_taken, sub_error);
            still.clear();
            for (size_t i = 0; i < again.size(); ++i)
                if (sub_taken[i] && sub_error[i] > kTreesSettled) still.push_back(again[i]);
            again.swap(still);
        }
        for (size_t index = 0; index < joints.size(); ++index) {
            Joint& joint = joints[index];
            // A joint in a tree is done, but for its one-sided rows (a limit, the cone), which the trees leave.
            // Those, and every row of the other joints, one ROW at a time, each measured from the poses as they
            // are right now. A joint's rows share its bodies, so correcting one changes the error of the others:
            // closing a hinge's point turns the body, which changes its angle, and turning it back to a limit
            // moves the point. Corrected all together from one measurement they overshoot against each other and
            // never settle. Each row is measured afresh, and only that row is built. (Not once more after the
            // last: a joint has the same rows whatever its pose, except that a ball joint's cone row is left out
            // while the two axes are exactly in line, and a row added that way would have no error to correct.)
            // The exception is the three rows of an anchor point, corrected together (correct_anchor, above), in
            // the place they have among the joint's rows.
            const bool in_tree = taken[index] != 0;
            if (in_tree && !((joint.type == JointType::Ball && joint.enable_cone) ||
                             ((joint.type == JointType::Hinge || joint.type == JointType::Slider) && joint.enable_limit)))
                continue;
            const int anchor_rows = in_tree ? -1 : joint.type == JointType::Ball || joint.type == JointType::Hinge ? 0 : joint.type == JointType::Weld ? 3 : -1;
            if (anchor_rows >= 0 && !is_live(bodies, joint)) continue;
            int rows = anchor_rows >= 0 ? anchor_rows + 3 + 1 : 1;  // (until a prepare says how many there are)
            for (int k = 0; k < rows; ++k) {
                if (k == anchor_rows) {
                    correct_anchor(bodies, joint);
                    k += 2;
                    continue;
                }
                solver.prepare(bodies, std::span<Joint>(&joint, 1), 1, settings, /*for_velocity=*/false, k);
                rows = solver.row_count();
                if (solver.rows().empty()) break;
                JointRow& r = solver.rows_mutable()[0];
                if (r.error == 0 || r.mass <= 0) continue;
                if (in_tree && r.lower == -kInf && r.upper == kInf) continue;  // (the tree's)
                const bool angular_only = r.lin_a.length_sq() == 0 && r.lin_b.length_sq() == 0;
                const Real limit = angular_only ? kMaxAngularCorrection : kMaxLinearCorrection;
                const Real C = std::clamp(r.error, -limit, limit);
                apply_displacement(r, -C * r.mass);
            }
        }
    }
}

}  // namespace phys3d
