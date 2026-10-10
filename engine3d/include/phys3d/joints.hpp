#pragma once
// Stage 10, step 8: joints in 3D.
//
// The method is the 2D engine's (see phys/joints.hpp, which explains it in full): every joint is a few
// SCALAR VELOCITY ROWS of one shared form,
//
//     J v  =  lin_a . v_a + ang_a . w_a + lin_b . v_b + ang_b . w_b ,
//
// solved by an impulse along the row, with the accumulated impulse clamped (unbounded for a rigid row,
// >= 0 for a limit, +-max for a motor) and kept between steps to warm start the next. Position error is
// NOT fed back as velocity; a separate pass after the bodies have moved shifts positions directly
// (solve_joint_positions), one joint at a time.
//
// What is new in 3D:
//   - the angular part of a row is a vector, and "how hard is it to turn the body this way" goes through
//     the inertia tensor:  J M^-1 J^T = |lin_a|^2/m_a + ang_a . (I_a^-1 ang_a) + (same for b);
//   - a point has three coordinates, so pinning two points together takes three rows (a ball and socket
//     joint), and that leaves all three rotations free. A hinge has to take two of them away again;
//   - a relative rotation is no longer one number. A hinge measures its angle about its own axis only.
#include "body.hpp"

#include <phys/solver_settings.hpp>

#include <span>
#include <vector>

namespace phys3d {

using phys::SolverSettings;

enum class JointType {
    Distance,  // keep two anchor points a fixed distance apart (a rod, or a spring when soft)
    Ball,      // pin two bodies together at one point; they may turn any way about it (a shoulder)
    Hinge,     // pin them together at a point AND keep them turning only about one axis (a door, a wheel);
               // optional angle limit and motor
};

struct Joint {
    // One slot per row a joint can produce, so a row keeps its accumulated impulse from step to step even
    // when its neighbours (a limit that is not currently reached) come and go.
    //   Distance: 0.   Ball: 0-2.   Hinge: 0-2 the point, 3-4 the axis, 5 the motor, 6-7 the limits.
    static constexpr int kMaxRows = 8;

    JointType type = JointType::Distance;
    int a = -1, b = -1;  // body indices; -1 means "the fixed world" (anchors are then world coordinates)
    Vec3 local_anchor_a, local_anchor_b;
    // Jointed bodies do not collide with each other unless this is set.
    bool collide_connected = false;

    // Distance: target separation of the anchors. frequency_hz > 0 makes it a spring-damper instead of
    // a rigid rod.
    Real length = 0;
    Real frequency_hz = 0;
    Real damping_ratio = 0;

    // Hinge. The axis is stored in each body's own frame; the joint keeps the two copies parallel. `ref`
    // is a direction perpendicular to the axis, again one copy per body: the hinge angle is how far b's
    // copy has turned from a's about the axis, in (-pi, pi], zero where the joint was made.
    Vec3 local_axis_a{1, 0, 0}, local_axis_b{1, 0, 0};
    Vec3 local_ref_a{0, 1, 0}, local_ref_b{0, 1, 0};
    bool enable_limit = false;  // keep the angle in [lower, upper]; both must lie inside (-pi, pi)
    Real lower = 0, upper = 0;
    bool enable_motor = false;  // drive the angular speed of b relative to a, about the axis
    Real motor_speed = 0;
    Real max_motor = 0;  // largest torque the motor may use

    // Solver state carried between steps.
    Real impulses[kMaxRows] = {};

    // Factories take world-space anchors and compute the local ones (and the rest state) from the bodies'
    // CURRENT poses, so a joint made from a valid scene starts with zero error. -1 attaches to the world.
    static Joint distance(const std::vector<Body>& bodies, int a, int b, Vec3 world_anchor_a, Vec3 world_anchor_b);
    static Joint ball(const std::vector<Body>& bodies, int a, int b, Vec3 world_anchor);
    static Joint hinge(const std::vector<Body>& bodies, int a, int b, Vec3 world_anchor, Vec3 world_axis);

    // World-space position of each anchor now (for drawing and for measuring joint error).
    Vec3 world_anchor_a(const std::vector<Body>& bodies) const;
    Vec3 world_anchor_b(const std::vector<Body>& bodies) const;
    // Hinge only: the axis as body a carries it, and the current angle.
    Vec3 world_axis(const std::vector<Body>& bodies) const;
    Real hinge_angle(const std::vector<Body>& bodies) const;
};

// One scalar constraint row, fully prepared for this step. (Public so tests can check the Jacobians.)
struct JointRow {
    Body* a = nullptr;
    Body* b = nullptr;
    Vec3 lin_a, ang_a, lin_b, ang_b;
    Vec3 inv_i_ang_a, inv_i_ang_b;  // I^-1 ang, world frame: the spin one unit of impulse gives each body
    Real bias = 0;
    Real error = 0;             // position error this row measures (0 if none); removed by solve_joint_positions
    Real softness = 0;          // "gamma": 0 for a rigid row
    Real mass = 0;              // 1 / (J M^-1 J^T + softness)
    Real lower = 0, upper = 0;  // bounds on the ACCUMULATED impulse
    Real impulse = 0;
    Real* store = nullptr;      // where to save it for the next step (the joint's slot)

    // J v for the bodies' current velocities.
    Real velocity() const { return dot(lin_a, a->vel) + dot(ang_a, a->w) + dot(lin_b, b->vel) + dot(ang_b, b->w); }
};

class JointSolver {
public:
    // Builds this step's rows from the bodies' current poses. `joints` and `bodies` must not be resized
    // until store() has been called. `for_velocity` false builds rows for the position pass: stored
    // impulses are left untouched.
    void prepare(std::vector<Body>& bodies, std::span<Joint> joints, Real dt, const SolverSettings& settings,
                 bool for_velocity = true);
    void warm_start();     // re-applies last step's impulses
    void solve_velocity(); // one sweep over all rows
    void store();          // saves the accumulated impulses back into the joints

    const std::vector<JointRow>& rows() const { return rows_; }
    std::vector<JointRow>& rows_mutable() { return rows_; }

private:
    std::vector<JointRow> rows_;
    Body world_;  // stands in for body index -1: immovable, at the origin, unrotated
};

// Position correction, run after the bodies have moved: shifts positions and orientations directly to
// remove whatever joint error is left (at most 20 cm / 8 degrees per correction). Velocities are not
// touched, so it cannot add energy. Rows are corrected one at a time, each measured from the poses the
// previous correction left behind (Gauss-Seidel), within a joint as well as between joints.
void solve_joint_positions(std::vector<Body>& bodies, std::vector<Joint>& joints, const SolverSettings& settings);

}  // namespace phys3d
