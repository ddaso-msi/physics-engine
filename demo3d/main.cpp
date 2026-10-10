// Wireframe demo for the 3D engine (phys3d). Everything is drawn as lines through a perspective camera:
// boxes as their 12 edges, spheres as three circles fixed in the sphere's own frame (so you can see it spin),
// capsules as two rings, four side lines and arcs over the ends,
// contacts as a yellow dot with a green normal, joints in orange (hinge axes in yellow). Lines further away are dimmer; nothing is hidden.
//
//   left-drag  orbit the camera          wheel  zoom
//   1 tower   2 pyramid   3 mixed rain   4 slope (friction and rolling)   5 tumble (zero gravity)
//   6 chains, a bridge and a spring (ball and distance joints)   7 hinges: a door, a flail, a motor
//   Space  shoot a ball from the camera   F  drop 10 more random bodies
//   W  warm starting on/off   - / =  fewer / more solver sweeps   G  gyroscopic term: implicit / explicit / off
//   N  narrow phase: per-shape routines / GJK and EPA
//   B  broad phase: tree / brute force   S  sleeping on/off (sleepers are drawn dim)   K  continuous collision on/off
//   C  show contacts   H  hide/show text   Tab  pause   Backspace or R  reload the scene   Esc  quit
//
// Headless capture (for checking a scene without a window):
//   SDL_VIDEODRIVER=dummy demo3d --scene 2 --steps 240 --capture out.bmp
#include "camera.hpp"

#include <SDL3/SDL.h>
#include <phys/timestep.hpp>
#include <phys3d/world.hpp>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

using namespace phys3d;
using demo3d::Camera;
using demo3d::ScreenPoint;

namespace {

constexpr int kWidth = 960, kHeight = 640;
constexpr size_t kMaxBodies = 1200;  // a debug build gets slow well before this
constexpr Real kDt = 1.0f / 120.0f;

enum class Scene { Tower = 1, Pyramid, Rain, Slope, Tumble, Chains, Machines };
constexpr int kSceneCount = 7;
const char* scene_name(Scene s) {
    switch (s) {
        case Scene::Tower: return "tower";
        case Scene::Pyramid: return "pyramid";
        case Scene::Rain: return "mixed rain";
        case Scene::Slope: return "slope";
        case Scene::Chains: return "chains and a bridge (ball joints, a spring)";
        case Scene::Machines: return "hinges: door, flail, motor";
        default: return "tumble (zero gravity)";
    }
}
const char* gyro_name(Gyroscopic g) {
    return g == Gyroscopic::Implicit ? "implicit" : g == Gyroscopic::Explicit ? "explicit" : "off";
}

struct Lcg {
    std::uint32_t s = 1;
    float next() {
        s = s * 1664525u + 1013904223u;
        return static_cast<float>(s >> 8) / static_cast<float>(1u << 24);
    }
    float range(float lo, float hi) { return lo + (hi - lo) * next(); }
};

Body crate(Vec3 pos, Vec3 half = {0.5f, 0.5f, 0.5f}, Quat q = {}) {
    Body b = Body::solid_box(half, 1, pos, q);
    b.restitution = 0.1f;
    return b;
}

void add_random_bodies(World& w, Lcg& rng, int count, float height) {
    for (int i = 0; i < count && w.bodies.size() < kMaxBodies; ++i) {
        const Vec3 pos{rng.range(-2.5f, 2.5f), height + static_cast<float>(i) * 0.9f, rng.range(-2.5f, 2.5f)};
        if (i % 3 == 0) {
            Body b = Body::solid_sphere(rng.range(0.3f, 0.55f), 1, pos);
            b.restitution = 0.3f;
            w.add(b);
        } else if (i % 3 == 1) {
            const Quat q = Quat::from_axis_angle({rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)}, rng.range(0, 3));
            Body b = Body::solid_capsule(rng.range(0.25f, 0.6f), rng.range(0.2f, 0.35f), 1, pos, q);
            b.restitution = 0.1f;
            w.add(b);
        } else {
            const Quat q = Quat::from_axis_angle({rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)}, rng.range(0, 3));
            w.add(crate(pos, {rng.range(0.25f, 0.6f), rng.range(0.25f, 0.6f), rng.range(0.25f, 0.6f)}, q));
        }
    }
}

void load_scene(Scene scene, World& w, Camera& cam, Lcg& rng) {
    w.truncate(0);
    w.gravity = {0, -9.81f, 0};
    rng.s = 1;
    cam = Camera{};
    Body floor = Body::fixed_box({20, 0.5f, 20}, {0, -0.5f, 0});  // top face at y = 0
    floor.restitution = 0;
    w.add(floor);

    switch (scene) {
        case Scene::Tower:
            for (int i = 0; i < 8; ++i) w.add(crate({0, 0.5f + static_cast<float>(i) * 1.001f, 0}));
            cam.target = {0, 3.5f, 0};
            cam.distance = 16;
            break;
        case Scene::Pyramid:
            for (int layer = 0; layer < 4; ++layer) {
                const int m = 4 - layer;
                for (int i = 0; i < m; ++i)
                    for (int k = 0; k < m; ++k)
                        w.add(crate({(static_cast<float>(i) - static_cast<float>(m - 1) * 0.5f) * 1.02f, 0.5f + static_cast<float>(layer) * 1.002f,
                                     (static_cast<float>(k) - static_cast<float>(m - 1) * 0.5f) * 1.02f}));
            }
            cam.target = {0, 1.5f, 0};
            cam.distance = 12;
            break;
        case Scene::Rain:
            add_random_bodies(w, rng, 40, 2.0f);
            cam.target = {0, 2, 0};
            cam.distance = 18;
            break;
        case Scene::Slope: {
            const Real theta = 0.35f;
            const Quat tilt = Quat::from_axis_angle({0, 0, 1}, theta);
            Body ramp = Body::fixed_box({6, 0.2f, 4}, {0, 2.5f, 0}, tilt);
            ramp.friction = 0.6f;
            ramp.restitution = 0;
            w.add(ramp);
            const Vec3 n = rotate(tilt, {0, 1, 0}), along = rotate(tilt, {1, 0, 0});
            const Vec3 top = Vec3{0, 2.5f, 0} + along * 4.5f + n * 0.2f;  // near the upper end of the ramp
            Body grippy = crate(top + n * 0.5f + Vec3{0, 0, -2.5f}, {0.5f, 0.5f, 0.5f}, tilt);
            grippy.friction = 0.9f;  // sqrt(0.9 x 0.6) = 0.73 > tan(theta): it stays
            w.add(grippy);
            Body icy = crate(top + n * 0.5f, {0.5f, 0.5f, 0.5f}, tilt);
            icy.friction = 0.02f;    // it slides
            w.add(icy);
            Body roller = Body::solid_sphere(0.5f, 1, top + n * 0.5f + Vec3{0, 0, 2.5f});
            roller.friction = 0.9f;  // it rolls
            roller.restitution = 0;
            w.add(roller);
            // A capsule lying across the slope rolls too, on two contact points.
            Body log = Body::solid_capsule(0.6f, 0.35f, 1, top + n * 0.35f + Vec3{-2.5f, 0, 0} - along * 2.5f + Vec3{0, 0, 0},
                                           tilt * Quat::from_axis_angle({1, 0, 0}, kPi / 2));
            log.friction = 0.9f;
            log.restitution = 0;
            w.add(log);
            cam.target = {0, 2, 0};
            cam.distance = 17;
            cam.yaw = 0.35f;
            break;
        }
        case Scene::Tumble: {
            // Three identical bricks, each spun about a different one of its axes with a 1% wobble. About the
            // longest and shortest axes the spin is steady; about the middle one it keeps flipping over.
            w.gravity = {};
            for (int axis = 0; axis < 3; ++axis) {
                Body b = Body::solid_box({1.0f, 0.6f, 0.2f}, 1, {(static_cast<float>(axis) - 1) * 3.2f, 3, 0});
                Vec3 spin{0.05f, 0.05f, 0.05f};
                spin[axis] = 5;
                b.w = spin;
                w.add(b);
            }
            cam.target = {0, 3, 0};
            cam.distance = 13;
            cam.pitch = 0.2f;
            break;
        }
        case Scene::Chains: {
            // A chain hung by one end with a heavy ball on the other, started out sideways so it swings.
            int previous = -1;
            const Vec3 top{-5, 7, 0};
            for (int i = 0; i < 8; ++i) {
                const int link = w.add(crate(top + Vec3{0.3f + 0.6f * static_cast<float>(i), 0, 0}, {0.3f, 0.06f, 0.06f}));
                w.add_joint(Joint::ball(w.bodies, previous, link, top + Vec3{0.6f * static_cast<float>(i), 0, 0}));
                previous = link;
            }
            const int weight = w.add(Body::solid_sphere(0.45f, 1, top + Vec3{4.8f + 0.45f, 0, 0}));
            w.add_joint(Joint::ball(w.bodies, previous, weight, top + Vec3{4.8f, 0, 0}));

            // A plank bridge pinned to the world at both ends, with a load dropped on it.
            previous = -1;
            const Vec3 left{-3.5f, 2.5f, 3};
            for (int i = 0; i < 7; ++i) {
                const int plank = w.add(crate(left + Vec3{0.5f + static_cast<float>(i), 0, 0}, {0.5f, 0.06f, 0.6f}));
                for (float side : {-0.5f, 0.5f})  // two pins per joint, so the planks cannot twist freely
                    w.add_joint(Joint::ball(w.bodies, previous, plank, left + Vec3{static_cast<float>(i), 0, side}));
                previous = plank;
            }
            for (float side : {-0.5f, 0.5f}) w.add_joint(Joint::ball(w.bodies, previous, -1, left + Vec3{7, 0, side}));
            w.add(crate(left + Vec3{3.5f, 1.5f, 0}, {0.4f, 0.4f, 0.4f}));

            // A ball on a spring.
            const int bob = w.add(Body::solid_sphere(0.4f, 1, {4, 4.5f, -2}));
            Joint spring = Joint::distance(w.bodies, -1, bob, {4, 7, -2}, {4, 4.5f, -2});
            spring.frequency_hz = 1;
            spring.damping_ratio = 0.05f;
            w.add_joint(spring);
            w.bodies[static_cast<size_t>(bob)].vel = {1.5f, 0, 0};

            cam.target = {0, 3.5f, 0};
            cam.distance = 19;
            cam.yaw = 0.3f;
            break;
        }
        case Scene::Machines: {
            // A door on a vertical hinge that opens 100 degrees either way, pushed by a ball.
            w.add(Body::fixed_box({0.1f, 1.5f, 0.1f}, {-5, 1.5f, 0}));
            const int door = w.add(crate({-3.9f, 1.5f, 0}, {1.0f, 1.3f, 0.06f}));
            Joint hinge = Joint::hinge(w.bodies, -1, door, {-4.9f, 1.5f, 0}, {0, 1, 0});
            hinge.enable_limit = true;
            hinge.lower = -1.75f;
            hinge.upper = 1.75f;
            w.add_joint(hinge);
            Body push = Body::solid_sphere(0.4f, 1, {-3.4f, 1.6f, 4});
            push.vel = {0, 1, -9};
            w.add(push);

            // A two-part flail: an arm on a horizontal hinge, a second arm hinged to it at right angles.
            const Vec3 pivot{0, 6, 0};
            const int upper = w.add(crate(pivot + Vec3{1, 0, 0}, {1.0f, 0.1f, 0.1f}));
            w.add_joint(Joint::hinge(w.bodies, -1, upper, pivot, {0, 0, 1}));
            const int lower = w.add(crate(pivot + Vec3{2, 0, 1}, {0.1f, 0.1f, 1.0f}));
            w.add_joint(Joint::hinge(w.bodies, upper, lower, pivot + Vec3{2, 0, 0}, {1, 0, 0}));

            // A motor turning a paddle through a pile of crates.
            w.add(Body::fixed_box({0.15f, 0.6f, 0.15f}, {5, 0.6f, 0}));
            const int paddle = w.add(crate({5, 1.4f, 0}, {2.0f, 0.15f, 0.15f}));
            Joint motor = Joint::hinge(w.bodies, -1, paddle, {5, 1.4f, 0}, {0, 1, 0});
            motor.enable_motor = true;
            motor.motor_speed = 1.5f;
            motor.max_motor = 60;
            w.add_joint(motor);
            for (int i = 0; i < 4; ++i) w.add(crate({5 + 1.4f * (i % 2 ? 1.0f : -1.0f), 0.7f, 1.2f * (i < 2 ? 1.0f : -1.0f)}, {0.3f, 0.7f, 0.3f}));

            cam.target = {0, 2.5f, 0};
            cam.distance = 19;
            cam.yaw = 0.5f;
            cam.pitch = 0.45f;
            break;
        }
    }
}

struct Painter {
    SDL_Renderer* ren;
    const Camera& cam;

    // A world-space line, dimmer the further away it is.
    void line(Vec3 a, Vec3 b, SDL_Color c) const {
        ScreenPoint pa, pb;
        Real depth = 0;
        if (!cam.project_segment(a, b, pa, pb, &depth)) return;
        const float fade = std::clamp(1.25f - depth / 45.0f, 0.25f, 1.0f);
        SDL_SetRenderDrawColor(ren, static_cast<Uint8>(static_cast<float>(c.r) * fade), static_cast<Uint8>(static_cast<float>(c.g) * fade),
                               static_cast<Uint8>(static_cast<float>(c.b) * fade), 255);
        SDL_RenderLine(ren, pa.x, pa.y, pb.x, pb.y);
    }

    void box(const Body& b, SDL_Color c) const {
        const Vec3 h = b.shape.half_extents;
        Vec3 corner[8];
        for (int i = 0; i < 8; ++i)
            corner[i] = apply(b.transform(), {(i & 1) ? h.x : -h.x, (i & 2) ? h.y : -h.y, (i & 4) ? h.z : -h.z});
        for (int i = 0; i < 8; ++i)
            for (int bit = 0; bit < 3; ++bit)
                if (!(i & (1 << bit))) line(corner[i], corner[i | (1 << bit)], c);  // each edge once
    }

    void sphere(const Body& b, SDL_Color c) const {
        constexpr int kSegments = 28;
        const Real r = b.shape.radius;
        for (int plane = 0; plane < 3; ++plane) {
            Vec3 prev;
            for (int i = 0; i <= kSegments; ++i) {
                const Real a = 2 * kPi * static_cast<Real>(i) / kSegments;
                Vec3 local;
                local[(plane + 1) % 3] = r * std::cos(a);
                local[(plane + 2) % 3] = r * std::sin(a);
                const Vec3 p = apply(b.transform(), local);
                if (i > 0) line(prev, p, c);
                prev = p;
            }
        }
    }

    // A capsule: a ring at each end of its axis, four lines along it, and two half circles over each end.
    void capsule(const Body& b, SDL_Color c) const {
        constexpr int kSegments = 20;
        const Real r = b.shape.radius, h = b.shape.half_length;
        for (const Real end : {-h, h}) {
            Vec3 prev;
            for (int i = 0; i <= kSegments; ++i) {
                const Real a = 2 * kPi * static_cast<Real>(i) / kSegments;
                const Vec3 p = apply(b.transform(), {r * std::cos(a), end, r * std::sin(a)});
                if (i > 0) line(prev, p, c);
                prev = p;
            }
            for (int plane = 0; plane < 2; ++plane) {  // the cap: half circles in the xy and zy planes
                for (int i = 0; i <= kSegments / 2; ++i) {
                    const Real a = kPi * static_cast<Real>(i) / (kSegments / 2);
                    const Real out = r * std::cos(a), up = (end < 0 ? -1 : 1) * r * std::sin(a);
                    const Vec3 p = apply(b.transform(), plane == 0 ? Vec3{out, end + up, 0} : Vec3{0, end + up, out});
                    if (i > 0) line(prev, p, c);
                    prev = p;
                }
            }
        }
        for (const Vec3& side : {Vec3{r, 0, 0}, Vec3{-r, 0, 0}, Vec3{0, 0, r}, Vec3{0, 0, -r}})
            line(apply(b.transform(), side + Vec3{0, -h, 0}), apply(b.transform(), side + Vec3{0, h, 0}), c);
    }

    void grid() const {
        const SDL_Color c{60, 64, 78, 255};
        for (int i = -10; i <= 10; i += 2) {
            line({static_cast<Real>(i), 0, -10}, {static_cast<Real>(i), 0, 10}, c);
            line({-10, 0, static_cast<Real>(i)}, {10, 0, static_cast<Real>(i)}, c);
        }
    }

    void dot(Vec3 p, SDL_Color c) const {
        const Vec3 v = cam.to_view(p);
        if (v.z < cam.near_plane) return;
        const ScreenPoint s = cam.to_screen(v);
        SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, 255);
        const SDL_FRect r{s.x - 2.5f, s.y - 2.5f, 5, 5};
        SDL_RenderFillRect(ren, &r);
    }
};

void draw_world(SDL_Renderer* ren, const World& w, const Camera& cam, bool show_contacts) {
    const Painter paint{ren, cam};
    paint.grid();
    for (size_t i = 1; i < w.bodies.size(); ++i) {  // body 0 is the floor, shown by the grid
        const Body& b = w.bodies[i];
        const SDL_Color c = b.type == BodyType::Static ? SDL_Color{150, 154, 170, 255}
                            : !b.awake                   ? SDL_Color{70, 100, 130, 255}  // asleep
                                                         : SDL_Color{120, 200, 255, 255};
        if (b.shape.type == Shape::Type::Box) paint.box(b, c);
        else if (b.shape.type == Shape::Type::Capsule) paint.capsule(b, c);
        else paint.sphere(b, c);
    }
    // Joints: a line from each body's centre to its anchor (the two anchors coincide unless the joint has
    // come apart, or is a rod or spring), and a short line along a hinge's axis.
    for (const Joint& j : w.joints) {
        const SDL_Color c{255, 140, 90, 255};
        const Vec3 pa = j.world_anchor_a(w.bodies), pb = j.world_anchor_b(w.bodies);
        if (j.type == JointType::Distance) {
            paint.line(pa, pb, c);
        } else {
            if (j.a >= 0) paint.line(w.bodies[static_cast<size_t>(j.a)].pos, pa, c);
            if (j.b >= 0) paint.line(w.bodies[static_cast<size_t>(j.b)].pos, pb, c);
        }
        if (j.type == JointType::Hinge) {
            const Vec3 axis = j.world_axis(w.bodies);
            paint.line(pa - axis * 0.4f, pa + axis * 0.4f, {255, 220, 120, 255});
        }
        paint.dot(pa, c);
    }
    if (show_contacts)
        for (const ContactPair& c : w.contacts())
            for (int k = 0; k < c.manifold.count; ++k) {
                const Vec3 p = c.manifold.points[k].point;
                paint.line(p, p + c.manifold.normal * 0.35f, {110, 220, 140, 255});
                paint.dot(p, {255, 210, 90, 255});
            }
}

}  // namespace

int main(int argc, char** argv) {
    // Optional scripted run: --scene N --steps N --capture file.bmp
    int arg_scene = 1, arg_steps = 0;
    const char* capture = nullptr;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!std::strcmp(argv[i], "--scene")) arg_scene = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--steps")) arg_steps = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--capture")) capture = argv[i + 1];
    }

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    if (!SDL_CreateWindowAndRenderer("physics-engine 3D: wireframe", kWidth, kHeight, 0, &window, &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    World world;
    Camera cam;
    Lcg rng;
    Scene scene = static_cast<Scene>(std::clamp(arg_scene, 1, kSceneCount));
    load_scene(scene, world, cam, rng);
    for (int i = 0; i < arg_steps; ++i) world.step(kDt);

    bool show_contacts = true, show_text = true, paused = false, dragging = false;
    phys::FixedTimestep timestep(kDt);
    Uint64 last = SDL_GetTicksNS();
    double step_ms = 0;
    bool running = true;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = false;
            if (e.type == SDL_EVENT_KEY_DOWN) {
                switch (e.key.key) {
                    case SDLK_ESCAPE: running = false; break;
                    case SDLK_1: case SDLK_2: case SDLK_3: case SDLK_4: case SDLK_5: case SDLK_6: case SDLK_7:
                        scene = static_cast<Scene>(static_cast<int>(e.key.key - SDLK_1) + 1);
                        load_scene(scene, world, cam, rng);
                        break;
                    case SDLK_R: case SDLK_BACKSPACE: load_scene(scene, world, cam, rng); break;
                    case SDLK_TAB: paused = !paused; break;
                    case SDLK_H: show_text = !show_text; break;
                    case SDLK_C: show_contacts = !show_contacts; break;
                    case SDLK_W: world.solver.warm_starting = !world.solver.warm_starting; break;
                    case SDLK_MINUS: world.solver.iterations = std::max(1, world.solver.iterations - 1); break;
                    case SDLK_EQUALS: world.solver.iterations = std::min(30, world.solver.iterations + 1); break;
                    case SDLK_G:
                        world.gyroscopic = world.gyroscopic == Gyroscopic::Implicit   ? Gyroscopic::Explicit
                                           : world.gyroscopic == Gyroscopic::Explicit ? Gyroscopic::Off
                                                                                       : Gyroscopic::Implicit;
                        break;
                    case SDLK_F: add_random_bodies(world, rng, 10, 6.0f); break;
                    case SDLK_S: world.allow_sleep = !world.allow_sleep; break;
                    case SDLK_K: world.continuous = !world.continuous; break;
                    case SDLK_B:
                        world.broadphase = world.broadphase == BroadphaseKind::DynamicTree ? BroadphaseKind::BruteForce : BroadphaseKind::DynamicTree;
                        break;
                    case SDLK_N:
                        world.narrowphase = world.narrowphase == NarrowPhase::Convex ? NarrowPhase::Specialised : NarrowPhase::Convex;
                        break;
                    case SDLK_SPACE:
                        if (world.bodies.size() < kMaxBodies) {
                            Body shot = Body::solid_sphere(0.35f, 3, cam.eye() + cam.forward() * 1.5f);
                            shot.vel = cam.forward() * 22;
                            shot.restitution = 0.2f;
                            world.add(shot);
                        }
                        break;
                    default: break;
                }
            }
            if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN && e.button.button == SDL_BUTTON_LEFT) dragging = true;
            if (e.type == SDL_EVENT_MOUSE_BUTTON_UP && e.button.button == SDL_BUTTON_LEFT) dragging = false;
            if (e.type == SDL_EVENT_MOUSE_MOTION && dragging) cam.orbit(-e.motion.xrel * 0.008f, e.motion.yrel * 0.008f);
            if (e.type == SDL_EVENT_MOUSE_WHEEL) cam.zoom(e.wheel.y > 0 ? 0.9f : 1.1f);
        }

        const Uint64 now = SDL_GetTicksNS();
        const Real frame = static_cast<Real>(static_cast<double>(now - last) * 1e-9);
        last = now;
        if (!paused && !capture) {
            const Uint64 t0 = SDL_GetTicksNS();
            const int steps = timestep.advance(frame, [&](Real dt) { world.step(dt); });
            if (steps > 0) step_ms = static_cast<double>(SDL_GetTicksNS() - t0) * 1e-6 / steps;
        }

        SDL_SetRenderDrawColor(renderer, 20, 22, 28, 255);
        SDL_RenderClear(renderer);
        draw_world(renderer, world, cam, show_contacts);
        if (show_text) {
            SDL_SetRenderDrawColor(renderer, 150, 154, 170, 255);
            SDL_RenderDebugText(renderer, 16, 12, "drag orbit | wheel zoom | 1-7 scene | Space shoot | F drop 10 | R reload | Tab pause | H hide text");
            SDL_RenderDebugTextFormat(renderer, 16, 28, "scene: %s | bodies %d | contact pairs %d | step %.2f ms%s", scene_name(scene),
                                      static_cast<int>(world.bodies.size()) - 1, static_cast<int>(world.contacts().size()), step_ms,
                                      paused ? " | PAUSED" : "");
            SDL_RenderDebugTextFormat(renderer, 16, 44, "W warm starting: %s | -/= sweeps: %d | G gyroscopic: %s | C contacts",
                                      world.solver.warm_starting ? "on" : "OFF", world.solver.iterations, gyro_name(world.gyroscopic));
            SDL_RenderDebugTextFormat(renderer, 16, 60, "B broad phase: %s | %llu box tests -> %d candidate pairs -> %d contacts | N narrow phase: %s",
                                      world.broadphase == BroadphaseKind::DynamicTree ? "tree" : "brute force",
                                      static_cast<unsigned long long>(world.stats().broadphase_tests),
                                      static_cast<int>(world.stats().candidate_pairs), static_cast<int>(world.stats().contacts),
                                      world.narrowphase == NarrowPhase::Convex ? "GJK/EPA" : "per shape");
            SDL_RenderDebugTextFormat(renderer, 16, 76, "S sleep: %s (%d awake, %d islands) | K continuous collision: %s (%d swept, %d stopped)",
                                      world.allow_sleep ? "on" : "OFF", static_cast<int>(world.stats().awake_bodies),
                                      static_cast<int>(world.stats().islands), world.continuous ? "on" : "OFF",
                                      static_cast<int>(world.stats().ccd_swept), static_cast<int>(world.stats().ccd_hits));
        }
        if (capture) {
            SDL_Surface* surf = SDL_RenderReadPixels(renderer, nullptr);
            const bool ok = surf && SDL_SaveBMP(surf, capture);
            SDL_DestroySurface(surf);
            if (!ok) SDL_Log("capture failed: %s", SDL_GetError());
            running = false;
        }
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
