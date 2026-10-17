#include <phys3d/solver.hpp>

#include "solving_order.hpp"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <unordered_map>

namespace phys3d {

void tangent_basis(Vec3 n, Vec3& t1, Vec3& t2) {
    // Cross n with whichever world axis it is least aligned with, so the result is never near zero.
    t1 = (std::fabs(n.x) > static_cast<Real>(0.57735) ? Vec3{n.y, -n.x, 0} : Vec3{0, n.z, -n.y}).normalized();
    t2 = cross(n, t1);
}

void transfer_impulses(const std::vector<ContactPair>& previous, std::vector<ContactPair>& current) {
    auto key = [](const ContactPair& p) {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.a)) << 32) | static_cast<std::uint32_t>(p.b);
    };
    std::unordered_map<std::uint64_t, const ContactPair*> by_pair;
    by_pair.reserve(previous.size());
    for (const ContactPair& p : previous) by_pair.emplace(key(p), &p);

    for (ContactPair& now : current) {
        const auto it = by_pair.find(key(now));
        if (it == by_pair.end()) continue;
        const Manifold& before = it->second->manifold;
        for (int k = 0; k < now.manifold.count; ++k) {
            ContactPoint& cp = now.manifold.points[k];
            cp.normal_impulse = cp.tangent_impulse[0] = cp.tangent_impulse[1] = 0;
            for (int j = 0; j < before.count; ++j) {
                if (before.points[j].id == cp.id) {
                    cp.normal_impulse = before.points[j].normal_impulse;
                    cp.tangent_impulse[0] = before.points[j].tangent_impulse[0];
                    cp.tangent_impulse[1] = before.points[j].tangent_impulse[1];
                    break;
                }
            }
        }
    }
}

namespace {

// One row of a contact: the impulse along `dir` (the normal or a tangent) at a point ra from a's centre of
// mass and rb from b's. Everything about the two bodies' turning that the sweeps need is worked out here,
// once, so that a sweep never multiplies by an inertia matrix or takes a cross product:
//   the speed along dir is   dot(vb - va, dir) + dot(wb, rb x dir) - dot(wa, ra x dir)
//   an impulse lambda along dir changes   va by -dir lambda / ma,  wa by -IA^-1 (ra x dir) lambda,
//                                         vb by +dir lambda / mb,  wb by +IB^-1 (rb x dir) lambda.
// (the same Jacobian form the joint rows use).
struct Row {
    Vec3 ja, jb;  // ra x dir, rb x dir
    Vec3 ka, kb;  // IA^-1 ja, IB^-1 jb
    Real mass = 0;  // 1 / (1/ma + 1/mb + ja . ka + jb . kb)

    void prepare(Vec3 ra, Vec3 rb, Vec3 dir, const Mat3& inv_ia, const Mat3& inv_ib, Real inv_ma, Real inv_mb) {
        ja = cross(ra, dir);
        jb = cross(rb, dir);
        ka = inv_ia * ja;
        kb = inv_ib * jb;
        const Real k = inv_ma + inv_mb + dot(ja, ka) + dot(jb, kb);
        mass = k > 0 ? 1 / k : 0;
    }
};

struct PointConstraint {
    Row normal, tangent[2];
    Real bias = 0;  // the normal speed we want to reach (>= 0)
    Real normal_impulse = 0;
    Real tangent_impulse[2] = {0, 0};
};

// The two bodies' velocities while one contact is being solved. They are read into these at the start and
// written back at the end, so the rows in between work in registers; and they are written back only to a
// body that can move. Every contact with the floor would otherwise write the floor's (unchanging) velocity
// and the next would have to wait to read it: a chain through memory from one end of the pile to the other.
struct Motion {
    Vec3 va, wa, vb, wb;
};

struct ContactConstraint {
    Body* a = nullptr;
    Body* b = nullptr;
    Real inv_ma = 0, inv_mb = 0;
    Vec3 normal, tangent[2];
    Real friction = 0;
    int count = 0;
    size_t source = 0;  // index of the ContactPair this came from (dormant pairs are skipped)
    PointConstraint points[Manifold::kMaxPoints];

    Motion load() const { return {a->vel, a->w, b->vel, b->w}; }
    void store(const Motion& m) const {
        if (a->type == BodyType::Dynamic) a->vel = m.va, a->w = m.wa;
        if (b->type == BodyType::Dynamic) b->vel = m.vb, b->w = m.wb;
    }
    static Real speed(const Motion& m, const Row& r, Vec3 dir) { return dot(m.vb - m.va, dir) + dot(m.wb, r.jb) - dot(m.wa, r.ja); }
    void push(Motion& m, const Row& r, Vec3 dir, Real lambda) const {
        m.va -= dir * (inv_ma * lambda);
        m.wa -= r.ka * lambda;
        m.vb += dir * (inv_mb * lambda);
        m.wb += r.kb * lambda;
    }
};

}  // namespace

void solve_constraints(std::vector<Body>& bodies, std::vector<ContactPair>& contacts, std::vector<Joint>& joints, Real dt,
                       const SolverSettings& settings, JointSolver* reuse) {
    // Joint rows are prepared first but change no velocities, so the contacts below still measure their
    // arrival speeds from the unmodified velocities.
    std::optional<JointSolver> fresh;
    JointSolver& joint_solver = reuse ? *reuse : fresh.emplace();
    joint_solver.prepare(bodies, joints, dt, settings);

    std::vector<ContactConstraint> constraints;
    constraints.reserve(contacts.size());

    // Solved colour by colour (see solving_order.hpp).
    const std::vector<size_t> order = solving_order(bodies, contacts.size(), [&](size_t i, int& a, int& b) {
        a = contacts[i].a, b = contacts[i].b;
        return !contacts[i].dormant;
    });
    for (const size_t source : order) {
        const ContactPair& pair = contacts[source];
        ContactConstraint c;
        c.source = source;
        c.a = &bodies[static_cast<size_t>(pair.a)];
        c.b = &bodies[static_cast<size_t>(pair.b)];
        c.normal = pair.manifold.normal;
        tangent_basis(c.normal, c.tangent[0], c.tangent[1]);
        c.friction = std::sqrt(c.a->friction * c.b->friction);
        c.count = pair.manifold.count;
        const Real restitution = std::max(c.a->restitution, c.b->restitution);

        c.inv_ma = c.a->inv_mass;
        c.inv_mb = c.b->inv_mass;
        const Mat3 inv_ia = c.a->inv_inertia_world(), inv_ib = c.b->inv_inertia_world();
        const Motion m = c.load();
        for (int k = 0; k < c.count; ++k) {
            const ContactPoint& cp = pair.manifold.points[k];
            PointConstraint& p = c.points[k];
            const Vec3 ra = cp.point - c.a->pos, rb = cp.point - c.b->pos;
            p.normal.prepare(ra, rb, c.normal, inv_ia, inv_ib, c.inv_ma, c.inv_mb);
            for (int t = 0; t < 2; ++t) p.tangent[t].prepare(ra, rb, c.tangent[t], inv_ia, inv_ib, c.inv_ma, c.inv_mb);

            // Bounce (fixed up front from the arrival speed) or push-out of any overlap, as in 2D.
            const Real approach = c.speed(m, p.normal, c.normal);
            const Real bounce = approach < -settings.restitution_threshold ? -restitution * approach : 0;
            const Real push_out = settings.baumgarte / dt * std::max(cp.depth - settings.slop, Real(0));
            p.bias = std::max(bounce, push_out);
            if (settings.warm_starting) {
                p.normal_impulse = cp.normal_impulse;
                p.tangent_impulse[0] = cp.tangent_impulse[0];
                p.tangent_impulse[1] = cp.tangent_impulse[1];
            }
        }
        constraints.push_back(c);
    }

    // Warm start, as a separate pass after every constraint has measured its arrival speed.
    if (settings.warm_starting) {
        joint_solver.warm_start();
        for (ContactConstraint& c : constraints) {
            Motion m = c.load();
            for (int k = 0; k < c.count; ++k) {
                const PointConstraint& p = c.points[k];
                c.push(m, p.normal, c.normal, p.normal_impulse);
                for (int t = 0; t < 2; ++t) c.push(m, p.tangent[t], c.tangent[t], p.tangent_impulse[t]);
            }
            c.store(m);
        }
    }

    for (int sweep = 0; sweep < settings.iterations; ++sweep) {
        joint_solver.solve_velocity(/*exactly=*/sweep == settings.iterations - 1);
        for (ContactConstraint& c : constraints) {
            Motion m = c.load();
            // Friction first; the normal constraint gets the last word in each sweep.
            for (int k = 0; k < c.count; ++k) {
                PointConstraint& p = c.points[k];
                // Solve the two tangent rows, then limit the LENGTH of the accumulated friction impulse.
                Real total[2];
                for (int t = 0; t < 2; ++t) total[t] = p.tangent_impulse[t] - c.speed(m, p.tangent[t], c.tangent[t]) * p.tangent[t].mass;
                const Real limit = c.friction * p.normal_impulse;
                const Real length_sq = total[0] * total[0] + total[1] * total[1];
                if (length_sq > limit * limit) {
                    const Real scale = length_sq > 0 ? limit / std::sqrt(length_sq) : 0;
                    total[0] *= scale;
                    total[1] *= scale;
                }
                for (int t = 0; t < 2; ++t) {
                    c.push(m, p.tangent[t], c.tangent[t], total[t] - p.tangent_impulse[t]);
                    p.tangent_impulse[t] = total[t];
                }
            }
            for (int k = 0; k < c.count; ++k) {
                PointConstraint& p = c.points[k];
                const Real approach = c.speed(m, p.normal, c.normal);
                Real lambda = -(approach - p.bias) * p.normal.mass;
                const Real total = std::max(p.normal_impulse + lambda, Real(0));  // push, never pull
                lambda = total - p.normal_impulse;
                p.normal_impulse = total;
                c.push(m, p.normal, c.normal, lambda);
            }
            c.store(m);
        }
    }

    joint_solver.store();
    for (size_t i = 0; i < constraints.size(); ++i)
        for (int k = 0; k < constraints[i].count; ++k) {
            ContactPoint& cp = contacts[constraints[i].source].manifold.points[k];
            cp.normal_impulse = constraints[i].points[k].normal_impulse;
            cp.tangent_impulse[0] = constraints[i].points[k].tangent_impulse[0];
            cp.tangent_impulse[1] = constraints[i].points[k].tangent_impulse[1];
        }
}

}  // namespace phys3d
