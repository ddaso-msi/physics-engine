#include <phys/solver.hpp>

#include <algorithm>

namespace phys {

namespace {

struct PointConstraint {
    Vec2 ra, rb;                 // contact point relative to each body's centre of mass
    Real normal_mass = 0;        // 1 / (effective inverse mass along the normal)
    Real tangent_mass = 0;       // same along the tangent
    Real bias = 0;               // the normal speed we want to reach (>= 0)
    Real normal_impulse = 0;     // accumulated over all sweeps
    Real tangent_impulse = 0;
};

struct ContactConstraint {
    Body* a = nullptr;
    Body* b = nullptr;
    Vec2 normal, tangent;
    Real friction = 0;
    int count = 0;
    PointConstraint points[2];
};

// Velocity of b's contact point relative to a's.
Vec2 relative_velocity(const ContactConstraint& c, const PointConstraint& p) {
    return c.b->vel + cross(c.b->w, p.rb) - c.a->vel - cross(c.a->w, p.ra);
}

// Equal and opposite impulse P at the contact point: +P on b, -P on a.
void apply_impulse(ContactConstraint& c, const PointConstraint& p, Vec2 impulse) {
    c.a->vel -= impulse * c.a->inv_mass;
    c.a->w -= c.a->inv_inertia * cross(p.ra, impulse);
    c.b->vel += impulse * c.b->inv_mass;
    c.b->w += c.b->inv_inertia * cross(p.rb, impulse);
}

// How hard it is to accelerate the contact point along `dir`: the linear inverse masses plus the
// rotational part, where a push at lever arm r only spins the body by (r x dir)^2 / I.
Real effective_inverse_mass(const ContactConstraint& c, const PointConstraint& p, Vec2 dir) {
    const Real ra = cross(p.ra, dir), rb = cross(p.rb, dir);
    return c.a->inv_mass + c.b->inv_mass + ra * ra * c.a->inv_inertia + rb * rb * c.b->inv_inertia;
}

}  // namespace

void solve_contacts(std::vector<Body>& bodies, const std::vector<ContactPair>& contacts, Real dt,
                    const SolverSettings& settings) {
    std::vector<ContactConstraint> constraints;
    constraints.reserve(contacts.size());

    // Build the constraints from the velocities and overlaps we have right now.
    for (const ContactPair& pair : contacts) {
        ContactConstraint c;
        c.a = &bodies[static_cast<size_t>(pair.a)];
        c.b = &bodies[static_cast<size_t>(pair.b)];
        c.normal = pair.manifold.normal;
        c.tangent = {c.normal.y, -c.normal.x};
        c.friction = std::sqrt(c.a->friction * c.b->friction);  // either slippery surface wins out
        c.count = pair.manifold.count;
        const Real restitution = std::max(c.a->restitution, c.b->restitution);  // bouncier one wins

        for (int k = 0; k < c.count; ++k) {
            const ContactPoint& cp = pair.manifold.points[k];
            PointConstraint& p = c.points[k];
            p.ra = cp.point - c.a->pos;
            p.rb = cp.point - c.b->pos;

            const Real kn = effective_inverse_mass(c, p, c.normal);
            const Real kt = effective_inverse_mass(c, p, c.tangent);
            p.normal_mass = kn > 0 ? 1 / kn : 0;
            p.tangent_mass = kt > 0 ? 1 / kt : 0;

            // Restitution: reflect the approach speed, scaled by e. Fixed once, up front, from the
            // speed before any impulse; recomputing it each sweep would bounce repeatedly.
            const Real approach = dot(relative_velocity(c, p), c.normal);
            const Real bounce = approach < -settings.restitution_threshold ? -restitution * approach : 0;
            // Baumgarte: turn overlap beyond the slop into a separating speed that removes it over
            // roughly 1/baumgarte frames.
            const Real push_out = settings.baumgarte / dt * std::max(cp.depth - settings.slop, Real(0));
            p.bias = std::max(bounce, push_out);
        }
        constraints.push_back(c);
    }

    for (int sweep = 0; sweep < settings.iterations; ++sweep) {
        for (ContactConstraint& c : constraints) {
            // Friction first: it is bounded by the normal impulse, and the normal constraint is the
            // more important one so it gets the last word in each sweep.
            for (int k = 0; k < c.count; ++k) {
                PointConstraint& p = c.points[k];
                const Real slide = dot(relative_velocity(c, p), c.tangent);
                Real lambda = -slide * p.tangent_mass;
                const Real limit = c.friction * p.normal_impulse;  // Coulomb cone
                const Real total = std::clamp(p.tangent_impulse + lambda, -limit, limit);
                lambda = total - p.tangent_impulse;
                p.tangent_impulse = total;
                apply_impulse(c, p, c.tangent * lambda);
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
}

}  // namespace phys
