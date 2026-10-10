#include <phys3d/solver.hpp>

#include <algorithm>
#include <cstdint>
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

struct PointConstraint {
    Vec3 ra, rb;  // contact point relative to each body's centre of mass
    Real normal_mass = 0;
    Real tangent_mass[2] = {0, 0};
    Real bias = 0;  // the normal speed we want to reach (>= 0)
    Real normal_impulse = 0;
    Real tangent_impulse[2] = {0, 0};
};

struct ContactConstraint {
    Body* a = nullptr;
    Body* b = nullptr;
    Mat3 inv_ia, inv_ib;  // world-frame inverse inertia, fixed for the duration of the solve
    Vec3 normal, tangent[2];
    Real friction = 0;
    int count = 0;
    size_t source = 0;  // index of the ContactPair this came from (dormant pairs are skipped)
    PointConstraint points[Manifold::kMaxPoints];
};

Vec3 relative_velocity(const ContactConstraint& c, const PointConstraint& p) {
    return c.b->vel + cross(c.b->w, p.rb) - c.a->vel - cross(c.a->w, p.ra);
}

// Equal and opposite impulse at the contact point: +impulse on b, -impulse on a.
void apply_impulse(ContactConstraint& c, const PointConstraint& p, Vec3 impulse) {
    c.a->vel -= impulse * c.a->inv_mass;
    c.a->w -= c.inv_ia * cross(p.ra, impulse);
    c.b->vel += impulse * c.b->inv_mass;
    c.b->w += c.inv_ib * cross(p.rb, impulse);
}

Real effective_inverse_mass(const ContactConstraint& c, const PointConstraint& p, Vec3 dir) {
    return c.a->inv_mass + c.b->inv_mass + dot(dir, cross(c.inv_ia * cross(p.ra, dir), p.ra)) +
           dot(dir, cross(c.inv_ib * cross(p.rb, dir), p.rb));
}

}  // namespace

void solve_constraints(std::vector<Body>& bodies, std::vector<ContactPair>& contacts, std::vector<Joint>& joints, Real dt,
                       const SolverSettings& settings) {
    // Joint rows are prepared first but change no velocities, so the contacts below still measure their
    // arrival speeds from the unmodified velocities.
    JointSolver joint_solver;
    joint_solver.prepare(bodies, joints, dt, settings);

    std::vector<ContactConstraint> constraints;
    constraints.reserve(contacts.size());

    for (size_t source = 0; source < contacts.size(); ++source) {
        const ContactPair& pair = contacts[source];
        if (pair.dormant) continue;
        ContactConstraint c;
        c.source = source;
        c.a = &bodies[static_cast<size_t>(pair.a)];
        c.b = &bodies[static_cast<size_t>(pair.b)];
        c.inv_ia = c.a->inv_inertia_world();
        c.inv_ib = c.b->inv_inertia_world();
        c.normal = pair.manifold.normal;
        tangent_basis(c.normal, c.tangent[0], c.tangent[1]);
        c.friction = std::sqrt(c.a->friction * c.b->friction);
        c.count = pair.manifold.count;
        const Real restitution = std::max(c.a->restitution, c.b->restitution);

        for (int k = 0; k < c.count; ++k) {
            const ContactPoint& cp = pair.manifold.points[k];
            PointConstraint& p = c.points[k];
            p.ra = cp.point - c.a->pos;
            p.rb = cp.point - c.b->pos;

            const Real kn = effective_inverse_mass(c, p, c.normal);
            p.normal_mass = kn > 0 ? 1 / kn : 0;
            for (int t = 0; t < 2; ++t) {
                const Real kt = effective_inverse_mass(c, p, c.tangent[t]);
                p.tangent_mass[t] = kt > 0 ? 1 / kt : 0;
            }

            // Bounce (fixed up front from the arrival speed) or push-out of any overlap, as in 2D.
            const Real approach = dot(relative_velocity(c, p), c.normal);
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
        for (ContactConstraint& c : constraints)
            for (int k = 0; k < c.count; ++k) {
                const PointConstraint& p = c.points[k];
                apply_impulse(c, p, c.normal * p.normal_impulse + c.tangent[0] * p.tangent_impulse[0] + c.tangent[1] * p.tangent_impulse[1]);
            }
    }

    for (int sweep = 0; sweep < settings.iterations; ++sweep) {
        joint_solver.solve_velocity();
        for (ContactConstraint& c : constraints) {
            // Friction first; the normal constraint gets the last word in each sweep.
            for (int k = 0; k < c.count; ++k) {
                PointConstraint& p = c.points[k];
                const Vec3 rv = relative_velocity(c, p);
                // Solve the two tangent rows, then limit the LENGTH of the accumulated friction impulse.
                Real total[2];
                for (int t = 0; t < 2; ++t) total[t] = p.tangent_impulse[t] - dot(rv, c.tangent[t]) * p.tangent_mass[t];
                const Real limit = c.friction * p.normal_impulse;
                const Real length_sq = total[0] * total[0] + total[1] * total[1];
                if (length_sq > limit * limit) {
                    const Real scale = length_sq > 0 ? limit / std::sqrt(length_sq) : 0;
                    total[0] *= scale;
                    total[1] *= scale;
                }
                const Vec3 delta = c.tangent[0] * (total[0] - p.tangent_impulse[0]) + c.tangent[1] * (total[1] - p.tangent_impulse[1]);
                p.tangent_impulse[0] = total[0];
                p.tangent_impulse[1] = total[1];
                apply_impulse(c, p, delta);
            }
            for (int k = 0; k < c.count; ++k) {
                PointConstraint& p = c.points[k];
                const Real approach = dot(relative_velocity(c, p), c.normal);
                Real lambda = -(approach - p.bias) * p.normal_mass;
                const Real total = std::max(p.normal_impulse + lambda, Real(0));  // push, never pull
                lambda = total - p.normal_impulse;
                p.normal_impulse = total;
                apply_impulse(c, p, c.normal * lambda);
            }
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
