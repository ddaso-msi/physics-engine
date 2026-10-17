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

// ---- Four contacts at once (step 24) ---------------------------------------------------------------
//
// Contacts of one colour share no moving body, so they can be solved side by side. Four at a time here, one in
// each lane of a 4-wide vector: the same arithmetic as ContactConstraint, done on four of everything. (A vector
// type the compiler knows, GCC's and Clang's vector_size, which becomes NEON or SSE: written as plain loops over
// four lanes, the compiler unrolled them but judged vectorising them not worth it.) A lane without a contact, or
// a point slot a contact does not have, has no mass and so moves nothing.
constexpr int kLanes = 4;
using F4 = float __attribute__((vector_size(16)));

struct V4 {
    F4 x, y, z;
};
inline V4 operator-(V4 a, V4 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline F4 dot(V4 a, V4 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline void set(V4& v, int l, Vec3 u) { v.x[l] = u.x, v.y[l] = u.y, v.z[l] = u.z; }
inline Vec3 get(const V4& v, int l) { return {v.x[l], v.y[l], v.z[l]}; }

struct WideRow {
    V4 ja, jb, ka, kb;
    F4 mass;
};

struct WidePoint {
    WideRow normal, tangent[2];
    F4 bias, normal_impulse, tangent_impulse[2];
};

struct WideConstraint {
    WideConstraint() {}  // (no clearing: pack() sets what is read, a third of the kilobytes in a contact of one point)
    int source[kLanes];  // index of the contact, -1 for an empty lane
    Body* a[kLanes];
    Body* b[kLanes];
    F4 inv_ma, inv_mb, friction;
    V4 normal, tangent[2];
    int points;  // the most points any lane has
    WidePoint point[Manifold::kMaxPoints];
};

void set_row(WideRow& w, int l, const Row& r) {
    set(w.ja, l, r.ja), set(w.jb, l, r.jb), set(w.ka, l, r.ka), set(w.kb, l, r.kb);
    w.mass[l] = r.mass;
}

// Packs up to four scalar constraints (lanes beyond `count` empty, standing on `nothing`, a fixed body). Only the
// point slots in use are cleared: empty lanes and missing points must have no mass, and finite numbers
// everywhere, for 0 times anything else to be 0.
void pack(WideConstraint& w, const ContactConstraint* const* c, int count, Body* nothing) {
    w.points = 0;
    for (int l = 0; l < count; ++l) w.points = std::max(w.points, c[l]->count);
    const F4 zero = {0, 0, 0, 0};
    const V4 zero3 = {zero, zero, zero};
    w.inv_ma = w.inv_mb = w.friction = zero;
    w.normal = w.tangent[0] = w.tangent[1] = zero3;
    for (int i = 0; i < w.points; ++i) {
        WidePoint& q = w.point[i];
        for (WideRow* r : {&q.normal, &q.tangent[0], &q.tangent[1]}) r->ja = r->jb = r->ka = r->kb = zero3, r->mass = zero;
        q.bias = q.normal_impulse = q.tangent_impulse[0] = q.tangent_impulse[1] = zero;
    }
    for (int l = 0; l < kLanes; ++l) {
        w.source[l] = -1;
        w.a[l] = w.b[l] = nothing;
        if (l >= count) continue;
        const ContactConstraint& k = *c[l];
        w.a[l] = k.a, w.b[l] = k.b;
        w.inv_ma[l] = k.inv_ma, w.inv_mb[l] = k.inv_mb, w.friction[l] = k.friction;
        set(w.normal, l, k.normal), set(w.tangent[0], l, k.tangent[0]), set(w.tangent[1], l, k.tangent[1]);
        for (int i = 0; i < k.count; ++i) {
            const PointConstraint& p = k.points[i];
            WidePoint& q = w.point[i];
            set_row(q.normal, l, p.normal);
            set_row(q.tangent[0], l, p.tangent[0]);
            set_row(q.tangent[1], l, p.tangent[1]);
            q.bias[l] = p.bias;
            q.normal_impulse[l] = p.normal_impulse;
            q.tangent_impulse[0][l] = p.tangent_impulse[0];
            q.tangent_impulse[1][l] = p.tangent_impulse[1];
        }
    }
}

struct WideMotion {
    V4 va, wa, vb, wb;
    F4 speed(const WideRow& r, const V4& dir) const { return dot(vb - va, dir) + dot(wb, r.jb) - dot(wa, r.ja); }
    void push(const WideRow& r, const V4& dir, F4 inv_ma, F4 inv_mb, F4 lambda) {
        const F4 la = inv_ma * lambda, lb = inv_mb * lambda;
        va.x -= dir.x * la, va.y -= dir.y * la, va.z -= dir.z * la;
        wa.x -= r.ka.x * lambda, wa.y -= r.ka.y * lambda, wa.z -= r.ka.z * lambda;
        vb.x += dir.x * lb, vb.y += dir.y * lb, vb.z += dir.z * lb;
        wb.x += r.kb.x * lambda, wb.y += r.kb.y * lambda, wb.z += r.kb.z * lambda;
    }
};

// One sweep of four contacts, exactly as the scalar sweep does one: friction first, then the normals.
void solve_wide(WideConstraint& c) {
    WideMotion m;
    for (int l = 0; l < kLanes; ++l) set(m.va, l, c.a[l]->vel), set(m.wa, l, c.a[l]->w), set(m.vb, l, c.b[l]->vel), set(m.wb, l, c.b[l]->w);
    const F4 zero = {0, 0, 0, 0};
    for (int k = 0; k < c.points; ++k) {
        WidePoint& p = c.point[k];
        F4 total0 = p.tangent_impulse[0] - m.speed(p.tangent[0], c.tangent[0]) * p.tangent[0].mass;
        F4 total1 = p.tangent_impulse[1] - m.speed(p.tangent[1], c.tangent[1]) * p.tangent[1].mass;
        // Limit the LENGTH of the accumulated friction impulse: scale it back to the limit where it is longer.
        const F4 limit = c.friction * p.normal_impulse;
        const F4 length_sq = total0 * total0 + total1 * total1;
        F4 scale = {1, 1, 1, 1};
        for (int l = 0; l < kLanes; ++l)
            if (length_sq[l] > limit[l] * limit[l]) scale[l] = length_sq[l] > 0 ? limit[l] / std::sqrt(length_sq[l]) : 0.0f;
        total0 *= scale;
        total1 *= scale;
        m.push(p.tangent[0], c.tangent[0], c.inv_ma, c.inv_mb, total0 - p.tangent_impulse[0]);
        m.push(p.tangent[1], c.tangent[1], c.inv_ma, c.inv_mb, total1 - p.tangent_impulse[1]);
        p.tangent_impulse[0] = total0;
        p.tangent_impulse[1] = total1;
    }
    for (int k = 0; k < c.points; ++k) {
        WidePoint& p = c.point[k];
        F4 total = p.normal_impulse - (m.speed(p.normal, c.normal) - p.bias) * p.normal.mass;
        total = total > zero ? total : zero;  // push, never pull
        m.push(p.normal, c.normal, c.inv_ma, c.inv_mb, total - p.normal_impulse);
        p.normal_impulse = total;
    }
    for (int l = 0; l < kLanes; ++l) {
        if (c.a[l]->type == BodyType::Dynamic) c.a[l]->vel = get(m.va, l), c.a[l]->w = get(m.wa, l);
        if (c.b[l]->type == BodyType::Dynamic) c.b[l]->vel = get(m.vb, l), c.b[l]->w = get(m.wb, l);
    }
}

// The accumulated impulses back into the contacts the lanes came from, for next step's warm start.
void unpack(const WideConstraint& w, std::vector<ContactPair>& contacts) {
    for (int l = 0; l < kLanes; ++l) {
        if (w.source[l] < 0) continue;
        Manifold& m = contacts[static_cast<size_t>(w.source[l])].manifold;
        for (int i = 0; i < m.count; ++i) {
            m.points[i].normal_impulse = w.point[i].normal_impulse[l];
            m.points[i].tangent_impulse[0] = w.point[i].tangent_impulse[0][l];
            m.points[i].tangent_impulse[1] = w.point[i].tangent_impulse[1][l];
        }
    }
}

// Re-applies last step's impulses, four contacts at a time.
void warm_start_wide(WideConstraint& c) {
    WideMotion m;
    for (int l = 0; l < kLanes; ++l) set(m.va, l, c.a[l]->vel), set(m.wa, l, c.a[l]->w), set(m.vb, l, c.b[l]->vel), set(m.wb, l, c.b[l]->w);
    for (int k = 0; k < c.points; ++k) {
        const WidePoint& p = c.point[k];
        m.push(p.normal, c.normal, c.inv_ma, c.inv_mb, p.normal_impulse);
        for (int t = 0; t < 2; ++t) m.push(p.tangent[t], c.tangent[t], c.inv_ma, c.inv_mb, p.tangent_impulse[t]);
    }
    for (int l = 0; l < kLanes; ++l) {
        if (c.a[l]->type == BodyType::Dynamic) c.a[l]->vel = get(m.va, l), c.a[l]->w = get(m.wa, l);
        if (c.b[l]->type == BodyType::Dynamic) c.b[l]->vel = get(m.vb, l), c.b[l]->w = get(m.wb, l);
    }
}

}  // namespace

void solve_constraints(std::vector<Body>& bodies, std::vector<ContactPair>& contacts, std::vector<Joint>& joints, Real dt,
                       const SolverSettings& settings, JointSolver* reuse) {
    // Joint rows are prepared first but change no velocities, so the contacts below still measure their
    // arrival speeds from the unmodified velocities.
    std::optional<JointSolver> fresh;
    JointSolver& joint_solver = reuse ? *reuse : fresh.emplace();
    joint_solver.prepare(bodies, joints, dt, settings);

    // Solved colour by colour (see solving_order.hpp).
    std::vector<size_t> starts;
    std::vector<size_t> order = solving_order(
        bodies, contacts.size(),
        [&](size_t i, int& a, int& b) {
            a = contacts[i].a, b = contacts[i].b;
            return !contacts[i].dormant;
        },
        &starts);
    // Within a colour the order is free, so contacts with the same number of points are put together: four of
    // them are solved for as many points as the most any of them has.
    for (int k = 0; k + 1 < kColours; ++k)
        std::stable_sort(order.begin() + static_cast<std::ptrdiff_t>(starts[static_cast<size_t>(k)]),
                         order.begin() + static_cast<std::ptrdiff_t>(starts[static_cast<size_t>(k) + 1]),
                         [&](size_t x, size_t y) { return contacts[x].manifold.count > contacts[y].manifold.count; });

    auto prepare = [&](size_t source, ContactConstraint& c) {
        const ContactPair& pair = contacts[source];
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
            p.normal_impulse = settings.warm_starting ? cp.normal_impulse : 0;
            p.tangent_impulse[0] = settings.warm_starting ? cp.tangent_impulse[0] : 0;
            p.tangent_impulse[1] = settings.warm_starting ? cp.tangent_impulse[1] : 0;
        }
    };

    // Each colour in fours (the last few with empty lanes), prepared straight into its lanes; except the last
    // colour, which may repeat a body and so is solved one contact at a time.
    Body nothing = Body::fixed({});
    std::vector<WideConstraint> wide;
    wide.reserve(order.size() / kLanes + kColours);
    for (int k = 0; k + 1 < kColours; ++k)
        for (size_t i = starts[static_cast<size_t>(k)]; i < starts[static_cast<size_t>(k) + 1]; i += kLanes) {
            const int count = static_cast<int>(std::min<size_t>(kLanes, starts[static_cast<size_t>(k) + 1] - i));
            ContactConstraint four[kLanes];
            const ContactConstraint* lanes[kLanes];
            for (int l = 0; l < count; ++l) prepare(order[i + static_cast<size_t>(l)], four[l]), lanes[l] = &four[l];
            wide.emplace_back();
            pack(wide.back(), lanes, count, &nothing);
            for (int l = 0; l < count; ++l) wide.back().source[l] = static_cast<int>(four[l].source);
        }
    std::vector<ContactConstraint> shared(order.size() - starts[kColours - 1]);
    for (size_t i = 0; i < shared.size(); ++i) prepare(order[starts[kColours - 1] + i], shared[i]);

    // Warm start, as a separate pass after every constraint has measured its arrival speed.
    if (settings.warm_starting) {
        joint_solver.warm_start();
        for (WideConstraint& w : wide) warm_start_wide(w);
        for (ContactConstraint& c : shared) {
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
        joint_solver.solve_velocity(/*exactly=*/false);
        for (WideConstraint& w : wide) solve_wide(w);
        for (ContactConstraint& c : shared) {
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

    // The joints that form trees, exactly, after everything: the last word. (In the last sweep and before the
    // contacts, as it was, a hard landing's contact impulses came after it and opened the joints, by up to 2 cm on a
    // rope bridge with crates dropped on it, and closed them only in the position pass; solved last, 0.1 mm. What
    // it costs the contacts is a little of their closing speed, which the next step's push-out handles.)
    joint_solver.solve_exactly();
    joint_solver.store();
    for (const WideConstraint& w : wide) unpack(w, contacts);
    for (const ContactConstraint& c : shared)
        for (int k = 0; k < c.count; ++k) {
            ContactPoint& cp = contacts[c.source].manifold.points[k];
            cp.normal_impulse = c.points[k].normal_impulse;
            cp.tangent_impulse[0] = c.points[k].tangent_impulse[0];
            cp.tangent_impulse[1] = c.points[k].tangent_impulse[1];
        }
}

}  // namespace phys3d
