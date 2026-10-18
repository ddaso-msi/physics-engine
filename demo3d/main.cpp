// Wireframe demo for the 3D engine (phys3d). Everything is drawn as lines through a perspective camera:
// boxes as their 12 edges, spheres as three circles fixed in the sphere's own frame (so you can see it spin),
// capsules as two rings, four side lines and arcs over the ends, hulls as their edges,
// contacts as a yellow dot with a green normal, joints in orange (hinge axes in yellow). Lines further away are dimmer; nothing is hidden.
//
//   left-drag  orbit the camera          wheel  zoom
//   1 tower   2 pyramid   3 mixed rain   4 slope (friction and rolling)   5 tumble (zero gravity)
//   6 chains, a bridge and a spring (ball and distance joints)   7 hinges: a door, a flail, a motor
//   8 a rope bridge (hinges, held at both ends: a loop)   9 ragdolls down a flight of stairs
//   0 welded rods swinging into each other, and a welded beam holding a weight out from a post
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

enum class Scene { Tower = 1, Pyramid, Rain, Slope, Tumble, Chains, Machines, RopeBridge, Ragdolls, Welds };
constexpr int kSceneCount = 10;
const char* scene_name(Scene s) {
    switch (s) {
        case Scene::Tower: return "tower";
        case Scene::Pyramid: return "pyramid";
        case Scene::Rain: return "mixed rain";
        case Scene::Slope: return "slope";
        case Scene::Chains: return "chains and a bridge (ball joints, a spring)";
        case Scene::Machines: return "hinges: door, flail, motor";
        case Scene::RopeBridge: return "rope bridge (hinges, a loop through the world)";
        case Scene::Ragdolls: return "ragdolls (ball joints with cones, hinges with limits)";
        case Scene::Welds: return "welded rods, and a welded beam";
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
        } else if (i % 6 == 4) {
            // A lump: the hull of a dozen random points. Or, every other time, a prism with 3 to 8 sides.
            std::vector<Vec3> points;
            if (i % 12 == 4) {
                for (int k = 0; k < 12; ++k) points.push_back({rng.range(-0.5f, 0.5f), rng.range(-0.5f, 0.5f), rng.range(-0.5f, 0.5f)});
            } else {
                const int sides = 3 + static_cast<int>(rng.next() * 6);
                for (int k = 0; k < sides; ++k) {
                    const float a = 2 * kPi * static_cast<float>(k) / static_cast<float>(sides);
                    points.push_back({0.45f * std::cos(a), -0.3f, 0.45f * std::sin(a)});
                    points.push_back({0.45f * std::cos(a), 0.3f, 0.45f * std::sin(a)});
                }
            }
            const Quat q = Quat::from_axis_angle({rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)}, rng.range(0, 3));
            if (std::optional<Body> b = Body::solid_hull(points.data(), static_cast<int>(points.size()), 1, pos, q)) {
                b->restitution = 0.1f;
                w.add(*b);
            }
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

// A ragdoll standing at `at` (the middle of its torso): a torso; a head on a neck; arms and legs of two capsules
// each. Shoulders, hips and neck are ball joints kept inside a cone; elbows and knees are hinges that bend one way.
// Every joint is part of one tree, solved exactly, and damped (step 29), as a body's joints are: without it nothing
// stops a limb swinging or the head rolling on its neck, and a ragdoll at rest never falls asleep. (The parts are
// light, a few grams to a few hundred at density 1, so the damping is small; at these values the three come to rest
// and fall asleep about six seconds after they are pushed.)
void add_ragdoll(World& w, Vec3 at, Vec3 velocity) {
    const size_t first = w.bodies.size();
    const int torso = w.add(crate(at, {0.22f, 0.32f, 0.12f}));
    const int head = w.add(Body::solid_sphere(0.15f, 1, at + Vec3{0, 0.5f, 0}));
    Joint neck = Joint::ball(w.bodies, torso, head, at + Vec3{0, 0.35f, 0});
    neck.limit_cone(w.bodies, {0, 1, 0}, 0.6f);
    neck.damping = 0.02f;
    w.add_joint(neck);
    // Limbs: [shoulder or hip, the way it hangs, upper length, lower length, radius, cone, which way the middle joint bends]
    struct Limb {
        Vec3 root;
        Real upper, lower, radius, cone, bend;
    };
    const Limb limbs[4] = {{{0.3f, 0.26f, 0}, 0.13f, 0.12f, 0.06f, 1.9f, 1},
                           {{-0.3f, 0.26f, 0}, 0.13f, 0.12f, 0.06f, 1.9f, 1},
                           {{0.12f, -0.36f, 0}, 0.18f, 0.17f, 0.08f, 1.2f, -1},
                           {{-0.12f, -0.36f, 0}, 0.18f, 0.17f, 0.08f, 1.2f, -1}};
    for (const Limb& limb : limbs) {
        const Vec3 root = at + limb.root;
        const Real upper_length = 2 * (limb.upper + limb.radius), lower_length = 2 * (limb.lower + limb.radius);
        const int upper = w.add(Body::solid_capsule(limb.upper, limb.radius, 1, root - Vec3{0, upper_length / 2, 0}));
        Joint ball = Joint::ball(w.bodies, torso, upper, root);
        ball.limit_cone(w.bodies, {0, -1, 0}, limb.cone);
        ball.damping = 0.02f;
        w.add_joint(ball);
        const Vec3 middle = root - Vec3{0, upper_length, 0};
        const int lower = w.add(Body::solid_capsule(limb.lower, limb.radius, 1, middle - Vec3{0, lower_length / 2, 0}));
        Joint hinge = Joint::hinge(w.bodies, upper, lower, middle, {1, 0, 0});
        hinge.enable_limit = true;  // an elbow bends forward, a knee back, each by up to about 140 degrees
        hinge.lower = limb.bend > 0 ? -2.4f : -0.05f;
        hinge.upper = limb.bend > 0 ? 0.05f : 2.4f;
        hinge.damping = 0.01f;
        w.add_joint(hinge);
    }
    for (size_t i = first; i < w.bodies.size(); ++i) w.bodies[i].vel = velocity;
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

            // A plank bridge pinned to the world at both ends, with a load dropped on it. It hangs with some slack, three
            // planks sloping down, one level, three up: pulled straight between its pins, it could only sag under the
            // load by stretching, and its joints stood 35 mm open.
            previous = -1;
            const Real slope = 0.25f;
            Vec3 at{-3.5f, 2.5f, 3};
            for (int i = 0; i < 7; ++i) {
                const Real angle = i < 3 ? -slope : i > 3 ? slope : 0;
                const Vec3 along{std::cos(angle), std::sin(angle), 0};
                const int plank = w.add(crate(at + along * 0.5f, {0.5f, 0.06f, 0.6f}, Quat::from_axis_angle({0, 0, 1}, angle)));
                for (float side : {-0.5f, 0.5f})  // two pins per joint, so the planks cannot twist freely
                    w.add_joint(Joint::ball(w.bodies, previous, plank, at + Vec3{0, 0, side}));
                previous = plank;
                at += along;
            }
            for (float side : {-0.5f, 0.5f}) w.add_joint(Joint::ball(w.bodies, previous, -1, at + Vec3{0, 0, side}));
            w.add(crate({(at.x - 3.5f) / 2, 4, 3}, {0.4f, 0.4f, 0.4f}));

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
            // A door on a vertical hinge that opens 100 degrees either way, pushed by a ball. The hinge stands 7 cm out
            // from the post: on the post's face, as it was, the door (12 cm thick) could not turn either way without an
            // edge going into the post, and once the contacts held it out (step 33) it stayed shut.
            w.add(Body::fixed_box({0.1f, 1.5f, 0.1f}, {-5, 1.5f, 0}));
            const int door = w.add(crate({-3.83f, 1.5f, 0}, {1.0f, 1.3f, 0.06f}));
            Joint hinge = Joint::hinge(w.bodies, -1, door, {-4.83f, 1.5f, 0}, {0, 1, 0});
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
        case Scene::RopeBridge: {
            // Sixteen planks hinged end to end and hung from two posts with slack, started as a V so it swings into
            // its sag. Fixed at both ends it is a loop through the world: cut at one end and put back exactly
            // through the rest (step 27). Then a few crates are dropped across it.
            const int planks = 16;
            const Real length = 0.55f, tilt = 0.32f;
            const Real span = static_cast<Real>(planks) * length * std::cos(tilt);
            Vec3 at{-span / 2, 4, 0};
            for (const Real x : {at.x, at.x + span}) w.add(Body::fixed_box({0.15f, 2, 0.15f}, {x, 2, -0.75f})), w.add(Body::fixed_box({0.15f, 2, 0.15f}, {x, 2, 0.75f}));
            int previous = -1;
            for (int i = 0; i < planks; ++i) {
                const Real angle = i < planks / 2 ? -tilt : tilt;
                const Vec3 along{std::cos(angle), std::sin(angle), 0};
                const int plank = w.add(crate(at + along * (length / 2), {length / 2, 0.04f, 0.6f}, Quat::from_axis_angle({0, 0, 1}, angle)));
                Joint hinge = Joint::hinge(w.bodies, previous, plank, at, {0, 0, 1});
                hinge.damping = 0.012f;  // (a rope's own stiffness, which settles the swing slowly; see step 29)
                w.add_joint(hinge);
                previous = plank;
                at += along * length;
            }
            Joint last = Joint::hinge(w.bodies, previous, -1, at, {0, 0, 1});
            last.damping = 0.012f;
            w.add_joint(last);
            for (int i = 0; i < 5; ++i) w.add(crate({-2.4f + 1.2f * static_cast<float>(i), 5.5f + 0.8f * static_cast<float>(i), 0}, {0.3f, 0.3f, 0.3f}));
            cam.target = {0, 2.5f, 0};
            cam.distance = 13;
            cam.yaw = 0.25f;
            cam.pitch = 0.15f;
            break;
        }
        case Scene::Ragdolls: {
            // A flight of six steps, and three ragdolls dropped at the top of it.
            for (int i = 0; i < 6; ++i) {
                Body step = Body::fixed_box({0.5f, 0.25f * static_cast<float>(i + 1), 2}, {3 - static_cast<float>(i), 0.25f * static_cast<float>(i + 1), 0});
                step.friction = 0.6f;
                w.add(step);
            }
            // (each pushed off the top step, down the flight; the lowest with its feet, a metre and a half below its
            // middle, just clear of the step)
            for (int i = 0; i < 3; ++i)
                add_ragdoll(w, {-1.7f, 4.6f + 1.6f * static_cast<float>(i), -1.1f + 1.1f * static_cast<float>(i)}, {2.2f + 0.4f * static_cast<float>(i), 0, 0});
            cam.target = {0.5f, 1.5f, 0};
            cam.distance = 12;
            cam.yaw = 0.6f;
            cam.pitch = 0.3f;
            break;
        }
        case Scene::Welds: {
            // Five rods of six to ten welded links, each hung from a ball joint and started leaning the same way:
            // being of different lengths they swing at different rates, drift out of step and knock into one
            // another. Each rod is a tree, solved exactly, so it stays straight; Gauss-Seidel alone let a rod like
            // this bend at its welds.
            for (int r = 0; r < 5; ++r) {
                const Vec3 top{-4 + 2 * static_cast<float>(r), 7.5f, -1.5f};
                const Vec3 lean = Vec3{-0.3f, -1, 0.12f * static_cast<float>(r - 2)}.normalized();
                const Quat q = Quat::from_axis_angle(cross({0, -1, 0}, lean), std::acos(-lean.y));
                int previous = -1;
                for (int i = 0; i < 6 + r; ++i) {
                    const int link = w.add(crate(top + lean * (0.3f + 0.6f * static_cast<float>(i)), {0.07f, 0.3f, 0.07f}, q));
                    w.add_joint(previous < 0 ? Joint::ball(w.bodies, -1, link, top) : Joint::weld(w.bodies, previous, link, top + lean * (0.6f * static_cast<float>(i))));
                    previous = link;
                }
            }
            // A beam of ten welded links sticking out of a post, a heavy ball welded to its tip.
            const Vec3 root{-4, 4, 4};
            w.add(Body::fixed_box({0.2f, 2, 0.2f}, {root.x - 0.2f, 2, root.z}));
            int previous = -1;
            for (int i = 0; i < 10; ++i) {
                const int link = w.add(crate(root + Vec3{0.3f + 0.6f * static_cast<float>(i), 0, 0}, {0.3f, 0.1f, 0.1f}));
                w.add_joint(Joint::weld(w.bodies, previous, link, root + Vec3{0.6f * static_cast<float>(i), 0, 0}));
                previous = link;
            }
            const int ball = w.add(Body::solid_sphere(0.45f, 4, root + Vec3{6.45f, 0, 0}));
            w.add_joint(Joint::weld(w.bodies, previous, ball, root + Vec3{6, 0, 0}));
            cam.target = {0, 3.5f, 0.5f};
            cam.distance = 17;
            cam.yaw = 0.35f;
            cam.pitch = 0.2f;
            break;
        }
    }
}

// How far apart the two anchors of the worst joint are, among those whose anchors should meet (ball, hinge, weld).
Real worst_joint_gap(const World& w) {
    Real worst = 0;
    for (const Joint& j : w.joints)
        if (j.type == JointType::Ball || j.type == JointType::Hinge || j.type == JointType::Weld)
            worst = std::max(worst, distance(j.world_anchor_a(w.bodies), j.world_anchor_b(w.bodies)));
    return worst;
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

    // A hull: every side of every face (each edge twice, once from each of its two faces).
    void hull(const Body& b, SDL_Color c) const {
        const Hull& h = *b.shape.hull;
        for (const Hull::Face& f : h.faces)
            for (int i = 0; i < f.count; ++i)
                if (h.corner(f, i) < h.corner(f, i + 1)) line(apply(b.transform(), h.point(f, i)), apply(b.transform(), h.point(f, i + 1)), c);
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
        else if (b.shape.type == Shape::Type::Hull) paint.hull(b, c);
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
                    case SDLK_1: case SDLK_2: case SDLK_3: case SDLK_4: case SDLK_5: case SDLK_6: case SDLK_7: case SDLK_8: case SDLK_9:
                        scene = static_cast<Scene>(static_cast<int>(e.key.key - SDLK_1) + 1);
                        load_scene(scene, world, cam, rng);
                        break;
                    case SDLK_0:
                        scene = Scene::Welds;
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
            if (e.type == SDL_EVENT_MOUSE_WHEEL) cam.scroll(e.wheel.y);
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
            SDL_RenderDebugText(renderer, 16, 12, "drag orbit | wheel zoom | 1-9, 0 scene | Space shoot | F drop 10 | R reload | Tab pause | H hide text");
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
            if (!world.joints.empty())
                SDL_RenderDebugTextFormat(renderer, 16, 92, "joints %d | worst gap between a joint's anchors: %.4f mm", static_cast<int>(world.joints.size()),
                                          static_cast<double>(worst_joint_gap(world)) * 1000);
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
