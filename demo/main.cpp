// Stage 5 demo: collision response sandbox.
//
// Click to drop a shape at the cursor. Things to try:
//   - drop a rubber ball and a wooden one: bounce height scales with restitution (e^2 per bounce)
//   - drop a crate on a ramp, then switch to ice: friction decides if it slides or holds
//   - press P for a 28-crate pyramid, then W to turn warm starting off and watch it slump, or
//     lower the solver sweeps with - and = while warm starting is on and see how little it needs
//
//   1 box   2 circle   3 hexagon   4 triangle     M cycle material (rubber / wood / ice)
//   W warm starting on/off    - / =  fewer / more solver sweeps     C show contacts
//   P pyramid level    R default level    X clear dynamic bodies    Esc quit
//
// Gallery (stage 9):
//   L  cycle through scenes: sandbox, pyramid, joints, Newton's cradle, rag doll, car, shooting range
//   S  sleeping on/off (sleeping bodies are drawn dimmed)     K  continuous collision on/off
//   cradle: Space lifts the first ball again           rag doll: click drops a rag doll
//   car: left/right arrows drive (release to brake)   range: click or Space fires a bullet at 70 m/s
//
// Recording:
//   H  hide/show all on-screen text     Backspace  reload the current scene     Tab  pause/resume
//
// Joints (stage 8):
//   J  joint level: rope bridge (hinges), rope pendulum and a crate on a spring (distance joints),
//      a motorised paddle (hinge + motor), a crate on a tilted rail with end stops (prismatic)
//   right mouse button: grab any dynamic body and drag it (a mouse joint: a capped soft spring)
//
// Broad phase (stage 7):
//   F  rain 100 random shapes in from the top, one at a time into free space; stops when the room is
//      full (watch the box-test count)
//   B  cycle broad phase: tree / sweep and prune / brute force
//   T  draw boxes: the tree's padded boxes (tree mode) or each body's tight box (other modes)
//   G  ghost mode: new shapes ignore other ghosts (collision filtering) but still hit everything else
#include <SDL3/SDL.h>
#include <phys/timestep.hpp>
#include <phys/world.hpp>
#include <scenes.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace phys;

namespace {

constexpr int kWidth = 960, kHeight = 640;
constexpr float kPixelsPerMeter = 60.0f;
constexpr float kWorldW = static_cast<float>(kWidth) / kPixelsPerMeter;
constexpr float kWorldH = static_cast<float>(kHeight) / kPixelsPerMeter;
constexpr size_t kMaxBodies = 1500;  // the solver and broad phase handle this; a debug build gets slow

SDL_FPoint to_screen(Vec2 w) { return {w.x * kPixelsPerMeter, (kWorldH - w.y) * kPixelsPerMeter}; }
Vec2 to_world(float sx, float sy) { return {sx / kPixelsPerMeter, kWorldH - sy / kPixelsPerMeter}; }

struct Material {
    const char* name;
    Real restitution, friction;
};
constexpr Material kMaterials[] = {{"rubber", 0.8f, 0.6f}, {"wood", 0.2f, 0.5f}, {"ice", 0.0f, 0.03f}};

enum class Kind { Box, Circle, Hexagon, Triangle };
const char* kind_name(Kind k) {
    switch (k) {
        case Kind::Box: return "box";
        case Kind::Circle: return "circle";
        case Kind::Hexagon: return "hexagon";
        default: return "triangle";
    }
}

Body make_body(Kind kind, Vec2 pos, Real angle, const Material& m) {
    Shape shape;
    switch (kind) {
        case Kind::Box: shape = Shape::make_polygon(Polygon::box(0.5f, 0.5f)); break;
        case Kind::Circle: shape = Shape::make_circle(0.5f); break;
        case Kind::Hexagon: {
            std::vector<Vec2> pts;
            for (int k = 0; k < 6; ++k) {
                float a = 2.0f * kPi * static_cast<float>(k) / 6.0f;
                pts.push_back({0.6f * std::cos(a), 0.6f * std::sin(a)});
            }
            shape = Shape::make_polygon(*Polygon::from_points(pts));
            break;
        }
        case Kind::Triangle: {
            const Vec2 tri[3] = {{0, 0}, {1.3f, 0}, {0.4f, 1.1f}};
            shape = Shape::make_polygon(*Polygon::from_points(tri));
            break;
        }
    }
    Body b(shape, pos, angle);
    b.restitution = m.restitution;
    b.friction = m.friction;
    return b;
}

Body make_static_box(Vec2 half, Vec2 pos, Real angle) {
    Body b(Shape::make_polygon(Polygon::box(half.x, half.y)), pos, angle, BodyType::Static);
    b.restitution = 0.0f;
    b.friction = 0.5f;
    return b;
}

void build_level(World& world, size_t& static_count, bool pyramid) {
    world.truncate(0);
    world.add(make_static_box({kWorldW * 0.5f, 0.25f}, {kWorldW * 0.5f, 0.25f}, 0));            // floor
    world.add(make_static_box({0.25f, kWorldH * 0.5f}, {0.25f, kWorldH * 0.5f}, 0));            // left wall
    world.add(make_static_box({0.25f, kWorldH * 0.5f}, {kWorldW - 0.25f, kWorldH * 0.5f}, 0));  // right wall
    world.add(make_static_box({kWorldW * 0.5f, 0.25f}, {kWorldW * 0.5f, kWorldH + 0.25f}, 0));  // ceiling, just out of view

    if (pyramid) {
        static_count = world.bodies.size();
        constexpr int kBase = 7;  // 7 + 6 + ... + 1 = 28 crates
        for (int row = 0; row < kBase; ++row)
            for (int i = 0; i < kBase - row; ++i) {
                const float x = kWorldW * 0.5f + (static_cast<float>(i) - static_cast<float>(kBase - row - 1) * 0.5f) * 1.02f;
                world.add(make_body(Kind::Box, {x, 1.0f + static_cast<float>(row) * 1.002f}, 0, kMaterials[1]));
            }
        return;
    }

    world.add(make_static_box({2.8f, 0.15f}, {3.6f, 6.2f}, -0.3f));  // ramp, slopes right
    world.add(make_static_box({2.8f, 0.15f}, {11.6f, 4.4f}, 0.3f));  // ramp, slopes left
    static_count = world.bodies.size();

    for (int i = 0; i < 4; ++i)  // a tower of crates
        world.add(make_body(Kind::Box, {8.0f, 1.0f + static_cast<float>(i) * 1.001f}, 0, kMaterials[1]));
    world.add(make_body(Kind::Circle, {2.6f, 8.5f}, 0, kMaterials[0]));
}

// A level with one of each joint kind. All bodies except the floor and walls are dynamic.
void build_joint_level(World& world, size_t& static_count) {
    world.truncate(0);
    world.add(make_static_box({kWorldW * 0.5f, 0.25f}, {kWorldW * 0.5f, 0.25f}, 0));
    world.add(make_static_box({0.25f, kWorldH * 0.5f}, {0.25f, kWorldH * 0.5f}, 0));
    world.add(make_static_box({0.25f, kWorldH * 0.5f}, {kWorldW - 0.25f, kWorldH * 0.5f}, 0));
    world.add(make_static_box({kWorldW * 0.5f, 0.25f}, {kWorldW * 0.5f, kWorldH + 0.25f}, 0));
    static_count = world.bodies.size();

    // 1. Rope bridge: 12 planks joined end to end by hinges, both ends pinned to the world. The planks
    //    are laid out along a sagging curve so every hinge starts exactly where it should.
    {
        constexpr int kLinks = 12;
        const Vec2 left{1.0f, 6.4f}, right{7.4f, 6.4f};
        std::vector<Vec2> pts;
        for (int i = 0; i <= kLinks; ++i) {
            const float t = static_cast<float>(i) / kLinks;
            pts.push_back({left.x + (right.x - left.x) * t, left.y - 1.4f * 4.0f * t * (1.0f - t)});
        }
        std::vector<int> link;
        for (int i = 0; i < kLinks; ++i) {
            const Vec2 d = pts[static_cast<size_t>(i) + 1] - pts[static_cast<size_t>(i)];
            Body plank(Shape::make_polygon(Polygon::box(d.length() * 0.5f, 0.07f)),
                       (pts[static_cast<size_t>(i)] + pts[static_cast<size_t>(i) + 1]) * 0.5f, std::atan2(d.y, d.x));
            plank.friction = 0.8f;
            plank.restitution = 0;
            link.push_back(world.add(plank));
        }
        world.add_joint(Joint::revolute(world.bodies, -1, link.front(), pts.front()));
        for (int i = 1; i < kLinks; ++i)
            world.add_joint(Joint::revolute(world.bodies, link[static_cast<size_t>(i) - 1], link[static_cast<size_t>(i)],
                                            pts[static_cast<size_t>(i)]));
        world.add_joint(Joint::revolute(world.bodies, link.back(), -1, pts.back()));
    }

    // 2. A rope pendulum: a ball on a rigid distance joint to a point in the world.
    {
        const Vec2 pivot{10.2f, 9.9f};
        const int ball = world.add(make_body(Kind::Circle, pivot + Vec2{3.0f * std::sin(1.0f), -3.0f * std::cos(1.0f)}, 0, kMaterials[0]));
        world.add_joint(Joint::distance(world.bodies, -1, ball, pivot, world.bodies[static_cast<size_t>(ball)].pos));
    }

    // 3. A crate on a spring: a soft distance joint, released from a stretch.
    {
        const Vec2 anchor{13.6f, 9.9f};
        const int crate = world.add(make_body(Kind::Box, {13.6f, 7.4f}, 0, kMaterials[1]));
        Joint spring = Joint::distance(world.bodies, -1, crate, anchor, world.bodies[static_cast<size_t>(crate)].pos);
        spring.length = 1.8f;
        spring.frequency_hz = 1.2f;
        spring.damping_ratio = 0.05f;
        world.add_joint(spring);
    }

    // 4. A motorised paddle: a bar on a hinge driven at 1.2 rad/s with a capped torque.
    {
        const Vec2 hub{12.6f, 3.4f};
        Body bar(Shape::make_polygon(Polygon::box(1.5f, 0.12f)), hub, 0.0f);
        bar.friction = 0.8f;
        const int idx = world.add(bar);
        Joint motor = Joint::revolute(world.bodies, -1, idx, hub);
        motor.enable_motor = true;
        motor.motor_speed = 1.2f;
        motor.max_motor = 800.0f;
        world.add_joint(motor);
    }

    // 5. A crate on a tilted rail with end stops: a prismatic joint with limits.
    {
        const Vec2 start{1.8f, 4.6f}, axis{1.0f, -0.3f};
        const int crate = world.add(make_body(Kind::Box, start, 0, kMaterials[1]));
        Joint rail = Joint::prismatic(world.bodies, -1, crate, start, axis);
        rail.enable_limit = true;
        rail.lower = 0.0f;
        rail.upper = 4.5f;
        world.add_joint(rail);
    }
}

void draw_joint(SDL_Renderer* ren, const World& world, const Joint& j) {
    const Vec2 pa = j.world_anchor_a(world.bodies), pb = j.world_anchor_b(world.bodies);
    const SDL_FPoint a = to_screen(pa), b = to_screen(pb);
    auto dot_at = [&](SDL_FPoint p, float r) {
        const SDL_FRect box{p.x - r, p.y - r, 2 * r, 2 * r};
        SDL_RenderFillRect(ren, &box);
    };
    switch (j.type) {
        case JointType::Distance:
            if (j.frequency_hz > 0) SDL_SetRenderDrawColor(ren, 200, 130, 255, 255);  // spring: violet
            else SDL_SetRenderDrawColor(ren, 255, 170, 80, 255);                      // rod: orange
            SDL_RenderLine(ren, a.x, a.y, b.x, b.y);
            dot_at(a, 3);
            dot_at(b, 3);
            break;
        case JointType::Revolute:
            SDL_SetRenderDrawColor(ren, 90, 220, 230, 255);
            dot_at(a, 4);
            break;
        case JointType::Prismatic: {
            const Vec2 axis = j.a < 0 ? j.local_axis_a : rotate(world.bodies[static_cast<size_t>(j.a)].q, j.local_axis_a);
            const Real lo = j.enable_limit ? j.lower : -3, hi = j.enable_limit ? j.upper : 3;
            const SDL_FPoint from = to_screen(pa + axis * lo), to = to_screen(pa + axis * hi);
            SDL_SetRenderDrawColor(ren, 110, 200, 120, 255);
            SDL_RenderLine(ren, from.x, from.y, to.x, to.y);
            dot_at(from, 3);
            dot_at(to, 3);
            break;
        }
        case JointType::Mouse: {
            const SDL_FPoint t = to_screen(j.target);
            SDL_SetRenderDrawColor(ren, 255, 90, 90, 255);
            SDL_RenderLine(ren, a.x, a.y, t.x, t.y);
            dot_at(t, 4);
            break;
        }
    }
}

enum class Scene { Sandbox, Pyramid, Joints, Cradle, Ragdoll, Car, Range };
constexpr int kSceneCount = 7;
const char* scene_name(Scene s) {
    switch (s) {
        case Scene::Sandbox: return "sandbox";
        case Scene::Pyramid: return "pyramid";
        case Scene::Joints: return "joints";
        case Scene::Cradle: return "Newton's cradle";
        case Scene::Ragdoll: return "rag doll";
        case Scene::Car: return "car";
        default: return "shooting range";
    }
}

// Floor, side walls and a ceiling just out of view, like the other levels.
void add_room(World& world) {
    world.add(make_static_box({kWorldW * 0.5f, 0.25f}, {kWorldW * 0.5f, 0.25f}, 0));
    world.add(make_static_box({0.25f, kWorldH * 0.5f}, {0.25f, kWorldH * 0.5f}, 0));
    world.add(make_static_box({0.25f, kWorldH * 0.5f}, {kWorldW - 0.25f, kWorldH * 0.5f}, 0));
    world.add(make_static_box({kWorldW * 0.5f, 0.25f}, {kWorldW * 0.5f, kWorldH + 0.25f}, 0));
}

struct Gallery {
    scenes::Cradle cradle;
    scenes::Car car;
    int ragdolls = 0;
    int rain_pending = 0;  // shapes F has asked for that have not been dropped yet
    int rain_steps = 0;    // physics steps since the last drop
    int rain_count = 0;    // shapes dropped so far: picks the next slot and shape
    int rain_blocked = 0;  // drops in a row that found every slot occupied
};

void load_scene(Scene scene, World& world, size_t& static_count, Gallery& g) {
    g = Gallery{};
    switch (scene) {
        case Scene::Sandbox: build_level(world, static_count, false); break;
        case Scene::Pyramid: build_level(world, static_count, true); break;
        case Scene::Joints: build_joint_level(world, static_count); break;
        case Scene::Cradle:
            world.truncate(0);
            add_room(world);
            static_count = world.bodies.size();
            g.cradle = scenes::build_cradle(world, {kWorldW * 0.5f, 9.8f});
            scenes::lift_cradle_ball(world, g.cradle, 0, -0.75f);
            break;
        case Scene::Ragdoll:
            world.truncate(0);
            add_room(world);
            world.add(make_static_box({1.6f, 0.2f}, {5.0f, 2.2f}, -0.35f));  // a ramp to tumble down
            static_count = world.bodies.size();
            scenes::build_ragdoll(world, {3.8f, 8.0f}, 1);
            g.ragdolls = 1;
            break;
        case Scene::Car:
            world.truncate(0);
            add_room(world);
            world.add(make_static_box({2.0f, 0.15f}, {10.5f, 0.55f}, 0.22f));  // a ramp up
            world.add(make_static_box({0.8f, 0.4f}, {5.5f, 0.7f}, 0));            // a bump
            static_count = world.bodies.size();
            g.car = scenes::build_car(world, {2.0f, 1.2f}, 1);
            break;
        case Scene::Range:
            world.truncate(0);
            add_room(world);
            world.add(make_static_box({0.03f, 4.5f}, {9.5f, 4.75f}, 0));  // a 6 cm wall
            static_count = world.bodies.size();
            for (int i = 0; i < 4; ++i)
                world.add(make_body(Kind::Box, {11.5f, 1.0f + static_cast<float>(i) * 1.001f}, 0, kMaterials[1]));  // resting on the floor (top at y = 0.5)
            break;
    }
}

// Rain (F): one shape every few steps, dropped into one of nine slots along the top of the room. The
// slots are wide enough that neighbours cannot touch whatever their angle, and a slot is only used when
// nothing is in it, so a shape never starts out overlapping anything: started overlapping, the contact
// push-out would fling it away at tens of m/s. A drop that finds every slot occupied (by shapes still
// falling clear, or by the pile) waits for the next turn; after a second of that the room is full and
// the rain stops.
constexpr int kRainSlots = 9;
constexpr int kRainInterval = 7;   // steps between drops: about 17 shapes a second
constexpr int kRainPatience = 17;  // blocked drops in a row before giving up: about a second

void rain_step(World& world, Gallery& g, std::uint32_t& seed) {
    if (g.rain_pending <= 0 || ++g.rain_steps < kRainInterval) return;
    g.rain_steps = 0;
    auto random = [&] {
        seed = seed * 1664525u + 1013904223u;
        return static_cast<float>(seed >> 8) / static_cast<float>(1u << 24);
    };
    const float slot_w = (kWorldW - 1.0f) / kRainSlots;
    const float reach = 0.71f;  // no shape reaches further than this from its centre (a box's corner)
    const float y = kWorldH - 0.9f;
    for (int attempt = 0; attempt < kRainSlots; ++attempt) {
        const int slot = (g.rain_count * 4 + attempt) % kRainSlots;  // stride 4: neighbours in time are apart in space
        const float x = 0.5f + slot_w * (static_cast<float>(slot) + 0.5f) + (random() - 0.5f) * 0.2f;
        const AABB space{{x - reach, y - reach}, {x + reach, y + reach}};
        bool is_free = true;
        for (const Body& other : world.bodies)
            if (overlap(space, compute_aabb(other))) {
                is_free = false;
                break;
            }
        if (!is_free) continue;
        const Kind k = (g.rain_count % 3 == 0) ? Kind::Circle : (g.rain_count % 3 == 1) ? Kind::Box : Kind::Hexagon;
        Body b = make_body(k, {x, y}, random() * 6.0f, kMaterials[1]);
        b.restitution = 0.1f;
        b.vel = {0.0f, -3.0f};
        world.add(b);
        ++g.rain_count;
        --g.rain_pending;
        g.rain_blocked = 0;
        return;
    }
    if (++g.rain_blocked >= kRainPatience) g.rain_pending = 0;  // the pile has reached the top
}

void fire_bullet(World& world, float y) {
    Body b(Shape::make_circle(0.12f), {1.0f, y}, 0, BodyType::Dynamic, 3.0f);
    b.restitution = 0.3f;
    b.vel = {70.0f, 0.0f};
    world.add(b);
}

void draw_aabb(SDL_Renderer* ren, const AABB& box) {
    const SDL_FPoint lo = to_screen(box.lo), hi = to_screen(box.hi);  // y is flipped: lo.y is the bottom
    const SDL_FRect r{lo.x, hi.y, hi.x - lo.x, lo.y - hi.y};
    SDL_RenderRect(ren, &r);
}

void draw_body(SDL_Renderer* ren, const Body& b, SDL_Color color) {
    SDL_SetRenderDrawColor(ren, color.r, color.g, color.b, 255);
    if (b.shape.type == Shape::Type::Circle) {
        constexpr int kSegments = 32;
        SDL_FPoint pts[kSegments + 1];
        const float r = b.shape.circle.radius;
        for (int i = 0; i <= kSegments; ++i) {
            float a = 2.0f * kPi * static_cast<float>(i) / kSegments;
            pts[i] = to_screen(b.pos + Vec2{r * std::cos(a), r * std::sin(a)});
        }
        SDL_RenderLines(ren, pts, kSegments + 1);
        SDL_FPoint c = to_screen(b.pos), tip = to_screen(b.pos + b.q.x_axis() * r);  // spoke shows spin
        SDL_RenderLine(ren, c.x, c.y, tip.x, tip.y);
    } else {
        const Polygon& p = b.shape.polygon;
        SDL_FPoint pts[Polygon::kMaxVertices + 1];
        for (int i = 0; i <= p.count; ++i) pts[i] = to_screen(apply(b.transform(), p.vertices[i % p.count]));
        SDL_RenderLines(ren, pts, p.count + 1);
    }
}

}  // namespace

int main() {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    if (!SDL_CreateWindowAndRenderer("physics-engine: collision response", kWidth, kHeight, 0, &window,
                                     &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    World world;
    size_t static_count = 0;
    Scene scene = Scene::Sandbox;
    Gallery gallery;
    load_scene(scene, world, static_count, gallery);

    Kind kind = Kind::Box;
    size_t material = 1;
    bool show_contacts = true;
    bool show_boxes = false;
    bool ghost_mode = false;
    std::uint32_t rain_seed = 1;
    int spawned = 0;
    int mouse_joint = -1;
    float last_mouse_y = 5.0f;  // index into world.joints of the joint dragging a body, or -1
    bool show_text = true;
    bool paused = false;

    FixedTimestep timestep(1.0f / 120.0f);
    Uint64 last = SDL_GetTicksNS();
    bool running = true;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = false;
            if (e.type == SDL_EVENT_KEY_DOWN) {
                switch (e.key.key) {
                    case SDLK_ESCAPE: running = false; break;
                    case SDLK_1: kind = Kind::Box; break;
                    case SDLK_2: kind = Kind::Circle; break;
                    case SDLK_3: kind = Kind::Hexagon; break;
                    case SDLK_4: kind = Kind::Triangle; break;
                    case SDLK_M: material = (material + 1) % (sizeof(kMaterials) / sizeof(kMaterials[0])); break;
                    case SDLK_C: show_contacts = !show_contacts; break;
                    case SDLK_X: mouse_joint = -1; gallery.rain_pending = 0; world.truncate(static_count); break;
                    case SDLK_R: mouse_joint = -1; scene = Scene::Sandbox; load_scene(scene, world, static_count, gallery); break;
                    case SDLK_P: mouse_joint = -1; scene = Scene::Pyramid; load_scene(scene, world, static_count, gallery); break;
                    case SDLK_J: mouse_joint = -1; scene = Scene::Joints; load_scene(scene, world, static_count, gallery); break;
                    case SDLK_L:
                        mouse_joint = -1;
                        scene = static_cast<Scene>((static_cast<int>(scene) + 1) % kSceneCount);
                        load_scene(scene, world, static_count, gallery);
                        break;
                    case SDLK_H: show_text = !show_text; break;
                    case SDLK_TAB: paused = !paused; break;
                    case SDLK_BACKSPACE:  // the current scene again, from its initial state
                        mouse_joint = -1;
                        rain_seed = 1;
                        spawned = 0;
                        load_scene(scene, world, static_count, gallery);
                        break;
                    case SDLK_S: world.allow_sleep = !world.allow_sleep; break;
                    case SDLK_K: world.continuous = !world.continuous; break;
                    case SDLK_SPACE:
                        if (scene == Scene::Cradle) scenes::lift_cradle_ball(world, gallery.cradle, 0, -0.75f);
                        if (scene == Scene::Range) fire_bullet(world, last_mouse_y);
                        break;
                    case SDLK_W: world.solver.warm_starting = !world.solver.warm_starting; break;
                    case SDLK_MINUS: world.solver.iterations = std::max(1, world.solver.iterations - 1); break;
                    case SDLK_EQUALS: world.solver.iterations = std::min(20, world.solver.iterations + 1); break;
                    case SDLK_T: show_boxes = !show_boxes; break;
                    case SDLK_G: ghost_mode = !ghost_mode; break;
                    case SDLK_B:
                        world.broadphase = world.broadphase == BroadphaseKind::DynamicTree      ? BroadphaseKind::SweepAndPrune
                                           : world.broadphase == BroadphaseKind::SweepAndPrune ? BroadphaseKind::BruteForce
                                                                                                : BroadphaseKind::DynamicTree;
                        break;
                    case SDLK_F:
                        gallery.rain_pending = std::min(gallery.rain_pending + 100,
                                                        static_cast<int>(kMaxBodies) - static_cast<int>(world.bodies.size()));
                        break;
                    default: break;
                }
            }
            // Right button: drag a body with a mouse joint, created on press and removed on release.
            if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN && e.button.button == SDL_BUTTON_RIGHT && mouse_joint < 0) {
                const Vec2 p = to_world(e.button.x, e.button.y);
                for (int i = static_cast<int>(world.bodies.size()) - 1; i >= 0; --i) {
                    const Body& b = world.bodies[static_cast<size_t>(i)];
                    if (b.type == BodyType::Dynamic && b.contains(p)) {
                        mouse_joint = world.add_joint(Joint::mouse(world.bodies, i, p));
                        break;
                    }
                }
            }
            if (e.type == SDL_EVENT_MOUSE_MOTION) last_mouse_y = to_world(e.motion.x, e.motion.y).y;
            if (e.type == SDL_EVENT_MOUSE_MOTION && mouse_joint >= 0 && mouse_joint < static_cast<int>(world.joints.size()))
                world.joints[static_cast<size_t>(mouse_joint)].target = to_world(e.motion.x, e.motion.y);
            if (e.type == SDL_EVENT_MOUSE_BUTTON_UP && e.button.button == SDL_BUTTON_RIGHT && mouse_joint >= 0) {
                if (mouse_joint < static_cast<int>(world.joints.size())) world.remove_joint(mouse_joint);
                mouse_joint = -1;
            }
            if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN && e.button.button == SDL_BUTTON_LEFT && scene == Scene::Range) {
                fire_bullet(world, to_world(e.button.x, e.button.y).y);
            } else if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN && e.button.button == SDL_BUTTON_LEFT && scene == Scene::Ragdoll &&
                       world.bodies.size() < kMaxBodies) {
                scenes::build_ragdoll(world, to_world(e.button.x, e.button.y), 1 + (gallery.ragdolls++ % 13));
            } else if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN && e.button.button == SDL_BUTTON_LEFT &&
                       world.bodies.size() < kMaxBodies) {
                // A small angle offset per spawn so stacked drops don't land perfectly aligned.
                const Real tilt = 0.17f * static_cast<Real>(spawned++ % 5) - 0.3f;
                Body b = make_body(kind, to_world(e.button.x, e.button.y), tilt, kMaterials[material]);
                if (ghost_mode) {  // group 2, collides only with group 1: not with other ghosts
                    b.category = 0b10;
                    b.mask = 0b01;
                }
                world.add(b);
            }
        }

        if (scene == Scene::Car) {
            const bool* keys = SDL_GetKeyboardState(nullptr);
            const Real speed = keys[SDL_SCANCODE_RIGHT] ? 12.0f : keys[SDL_SCANCODE_LEFT] ? -12.0f : 0.0f;
            scenes::set_car_throttle(world, gallery.car, speed, speed != 0 ? 25.0f : 40.0f);  // release = brake
        }

        Uint64 now = SDL_GetTicksNS();
        Real frame = static_cast<Real>(static_cast<double>(now - last) * 1e-9);
        last = now;
        // Paused: bank no time, so resuming carries on from here rather than catching up.
        if (!paused)
            timestep.advance(frame, [&](Real dt) {
                rain_step(world, gallery, rain_seed);
                world.step(dt);
            });

        SDL_SetRenderDrawColor(renderer, 20, 22, 28, 255);
        SDL_RenderClear(renderer);
        for (const Body& b : world.bodies)
            draw_body(renderer, b,
                      b.type == BodyType::Static ? SDL_Color{110, 114, 130, 255}
                      : !b.awake                 ? SDL_Color{70, 100, 130, 255}    // asleep
                      : b.category == 0b10       ? SDL_Color{190, 150, 255, 255}  // ghost
                                                 : SDL_Color{120, 200, 255, 255});

        for (const Joint& j : world.joints) draw_joint(renderer, world, j);

        if (show_boxes) {
            if (world.broadphase == BroadphaseKind::DynamicTree) {
                // Leaves are the padded boxes bodies are stored under; interior nodes are the
                // enclosing boxes, brighter the higher up the tree they are.
                world.broadphase_tree().for_each_node([&](const AABB& box, int height, bool leaf) {
                    if (leaf) SDL_SetRenderDrawColor(renderer, 90, 90, 40, 255);
                    else SDL_SetRenderDrawColor(renderer, 40, static_cast<Uint8>(std::min(60 + 22 * height, 255)), 90, 255);
                    draw_aabb(renderer, box);
                });
            } else {
                SDL_SetRenderDrawColor(renderer, 90, 90, 40, 255);
                for (const Body& b : world.bodies) draw_aabb(renderer, compute_aabb(b));
            }
        }

        if (show_contacts) {
            for (const ContactPair& c : world.contacts()) {
                for (int k = 0; k < c.manifold.count; ++k) {
                    const SDL_FPoint p = to_screen(c.manifold.points[k].point);
                    const SDL_FPoint tip = to_screen(c.manifold.points[k].point + c.manifold.normal * 0.4f);
                    SDL_SetRenderDrawColor(renderer, 110, 220, 140, 255);
                    SDL_RenderLine(renderer, p.x, p.y, tip.x, tip.y);
                    SDL_SetRenderDrawColor(renderer, 255, 210, 90, 255);
                    SDL_FRect dot{p.x - 2.5f, p.y - 2.5f, 5.0f, 5.0f};
                    SDL_RenderFillRect(renderer, &dot);
                }
            }
        }

        if (show_text) {
            SDL_SetRenderDrawColor(renderer, 150, 154, 170, 255);
            SDL_RenderDebugText(renderer, 70.0f, 12.0f, "click drop | 1-4 shape | M material | W warm start | -/= sweeps | P pyramid | R reset | X clear");
            const Material& m = kMaterials[material];
            SDL_RenderDebugTextFormat(renderer, 70.0f, 28.0f, "%s, %s (e=%.2f, mu=%.2f) | bodies %d | contacts %d",
                                      kind_name(kind), m.name, static_cast<double>(m.restitution),
                                      static_cast<double>(m.friction), static_cast<int>(world.bodies.size() - static_count),
                                      static_cast<int>(world.contacts().size()));
            SDL_RenderDebugTextFormat(renderer, 70.0f, 44.0f, "warm starting: %s | solver sweeps: %d | ghost mode: %s",
                                      world.solver.warm_starting ? "on" : "OFF", world.solver.iterations,
                                      ghost_mode ? "on" : "off");
            const char* bp_name = world.broadphase == BroadphaseKind::DynamicTree      ? "tree"
                                  : world.broadphase == BroadphaseKind::SweepAndPrune ? "sweep and prune"
                                                                                       : "brute force";
            SDL_RenderDebugTextFormat(renderer, 70.0f, 60.0f, "broad phase: %s | %llu box tests -> %d candidate pairs -> %d contacts",
                                      bp_name, static_cast<unsigned long long>(world.stats().broadphase_tests),
                                      static_cast<int>(world.stats().candidate_pairs), static_cast<int>(world.stats().contacts));
            SDL_RenderDebugText(renderer, 70.0f, 76.0f, "F rain 100 | B switch broad phase | T show boxes | G ghost mode");
            SDL_RenderDebugTextFormat(renderer, 70.0f, 92.0f, "J joint level | right-drag grabs a body | joints: %d",
                                      static_cast<int>(world.joints.size()));
            SDL_RenderDebugTextFormat(renderer, 70.0f, 108.0f,
                                      "L scene: %s | S sleep: %s (%d awake, %d islands) | K CCD: %s (%d swept, %d stopped)",
                                      scene_name(scene), world.allow_sleep ? "on" : "OFF",
                                      static_cast<int>(world.stats().awake_bodies), static_cast<int>(world.stats().islands),
                                      world.continuous ? "on" : "OFF", static_cast<int>(world.stats().ccd_swept),
                                      static_cast<int>(world.stats().ccd_hits));
            SDL_RenderDebugTextFormat(renderer, 70.0f, 124.0f, "H hide text | Backspace reload scene | Tab pause: %s",
                                      paused ? "PAUSED" : "running");
        }
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
