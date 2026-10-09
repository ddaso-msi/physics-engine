#pragma once
// Stage 8: joints.
//
// A joint removes some of the relative freedom of two bodies. Every joint here is expressed as a few
// SCALAR VELOCITY ROWS in one shared form, so the solver does not care which kind of joint a row came
// from:
//
//     J v  =  lin_a . v_a + ang_a * w_a + lin_b . v_b + ang_b * w_b
//
// is how fast the constrained quantity C is changing (the row is the derivative of C with respect to
// the two bodies' velocities). The solver picks an impulse lambda along the row, so the bodies feel
// (J^T lambda), and chooses lambda so that  J v + bias = 0:
//
//     lambda = -(J v + bias) / (J M^-1 J^T)         where M^-1 holds inverse masses and inertias.
//
// `bias` is the extra velocity a row asks for (a motor's target speed, a limit's allowance); position
// ERROR is kept out of the velocity pass entirely and removed by a separate pass after the bodies have
// moved (solve_joint_positions), which shifts positions directly. Feeding error back as velocity
// (Baumgarte) pumps energy into a joint, and with warm starting the pumping compounds from frame to
// frame: a chain with a hanging weight tears itself apart. A row may be
//   - an equality (accumulated impulse unbounded):    a rigid joint,
//   - an inequality (accumulated impulse >= 0):       a limit,
//   - bounded by +-max (accumulated impulse clamped): a motor with limited force,
//   - soft (a spring-damper): the effective mass is reduced and a restoring term added.
// Rows are rebuilt from the bodies' current poses every step; only each row's accumulated impulse is
// kept between steps, so it can warm start the next one.
#include "body.hpp"
#include "solver_settings.hpp"

#include <span>
#include <vector>

namespace phys {

enum class JointType {
    Distance,   // keep two anchor points a fixed distance apart (a rod, or a spring when soft)
    Revolute,   // pin two bodies together at one point (a hinge); optional angle limit and motor
    Prismatic,  // let body b slide along an axis fixed in body a, with no relative rotation (a rail)
    Mouse,      // pull a point on a body toward a target with a soft spring (for dragging)
};

struct Joint {
    // One slot per row a joint can produce, so a row keeps its accumulated impulse from step to step
    // even when its neighbours (a limit that is not currently touched) come and go.
    static constexpr int kMaxRows = 5;

    JointType type = JointType::Distance;
    int a = -1, b = -1;  // body indices; -1 means "the fixed world" (anchors are then world coordinates)
    Vec2 local_anchor_a, local_anchor_b;
    // Jointed bodies do not collide with each other unless this is set (ragdoll limbs would otherwise
    // fight the joint that holds them together).
    bool collide_connected = false;

    // Distance: target separation of the anchors. Distance and Mouse: frequency_hz > 0 makes the
    // joint a spring-damper instead of rigid (0 = rigid; ignored for Revolute/Prismatic).
    Real length = 0;
    Real frequency_hz = 0;
    Real damping_ratio = 0;

    // Revolute / Prismatic.
    Real reference_angle = 0;   // angle(b) - angle(a) at which the joint is relaxed
    Vec2 local_axis_a{1, 0};    // Prismatic: slide direction in body a's frame (unit)
    bool enable_limit = false;  // keep the angle (Revolute) or slide distance (Prismatic) in [lower, upper]
    Real lower = 0, upper = 0;
    bool enable_motor = false;  // drive the relative angular (Revolute) or sliding (Prismatic) speed
    Real motor_speed = 0;
    Real max_motor = 0;         // largest torque (Revolute) or force (Prismatic) the motor may use

    // Mouse.
    Vec2 target;
    Real max_force = 0;

    // Solver state carried between steps (see the header comment).
    Real impulses[kMaxRows] = {};

    // Factories take world-space anchors and compute the local anchors and rest state from the
    // bodies' CURRENT poses, so a joint made from a valid scene starts with zero error.
    // Use -1 for a (or b) to attach to the world.
    static Joint distance(const std::vector<Body>& bodies, int a, int b, Vec2 world_anchor_a, Vec2 world_anchor_b);
    static Joint revolute(const std::vector<Body>& bodies, int a, int b, Vec2 world_anchor);
    static Joint prismatic(const std::vector<Body>& bodies, int a, int b, Vec2 world_anchor, Vec2 world_axis);
    // Grabs `world_point` on body `body`; move it by changing `target`. max_force defaults to a few
    // times the body's weight times 10, enough to drag it around but not to launch it.
    static Joint mouse(const std::vector<Body>& bodies, int body, Vec2 world_point);

    // World-space position of each anchor now (for drawing and for measuring joint error).
    Vec2 world_anchor_a(const std::vector<Body>& bodies) const;
    Vec2 world_anchor_b(const std::vector<Body>& bodies) const;
};

// One scalar constraint row, fully prepared for this step. (Public so tests can check the Jacobians.)
struct JointRow {
    Body* a = nullptr;
    Body* b = nullptr;
    Vec2 lin_a, lin_b;
    Real ang_a = 0, ang_b = 0;
    Real bias = 0;
    Real error = 0;              // position error this row measures (0 if none); removed by solve_joint_positions
    Real softness = 0;           // "gamma": 0 for a rigid row
    Real mass = 0;               // 1 / (J M^-1 J^T + softness)
    Real lower = 0, upper = 0;   // bounds on the ACCUMULATED impulse
    Real impulse = 0;            // accumulated impulse so far
    Real* store = nullptr;       // where to save it for the next step (the joint's slot)
};

class JointSolver {
public:
    // Builds this step's rows from the bodies' current poses. `joints` and `bodies` must not be
    // resized until store() has been called.
    // `for_velocity` false builds rows for the position pass: stored impulses are left untouched.
    void prepare(std::vector<Body>& bodies, std::span<Joint> joints, Real dt, const SolverSettings& settings,
                 bool for_velocity = true);
    // Re-applies last step's impulses (if settings.warm_starting).
    void warm_start();
    // One sweep over all rows.
    void solve_velocity();
    // Saves the accumulated impulses back into the joints.
    void store();

    const std::vector<JointRow>& rows() const { return rows_; }
    std::vector<JointRow>& rows_mutable() { return rows_; }

private:
    std::vector<JointRow> rows_;
    Body world_;  // stands in for body index -1: immovable, at the origin, unrotated
};

// Position correction. Run after the bodies have been moved: shifts the bodies' positions and angles
// directly to remove whatever joint error is left (each correction moves them by at most 20 cm / 8
// degrees). Velocities are not touched, so this cannot add energy. Joints are corrected ONE AT A TIME,
// each from the poses the previous correction left behind (Gauss-Seidel): in a chain every link is
// shared by two joints, and correcting them all from the same starting poses would apply both
// corrections to the shared link and overshoot.
void solve_joint_positions(std::vector<Body>& bodies, std::vector<Joint>& joints, const SolverSettings& settings);

}  // namespace phys
