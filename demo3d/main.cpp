// Wireframe demo for the 3D engine (phys3d). Everything is drawn as lines through a perspective camera:
// boxes as their 12 edges, spheres as three circles fixed in the sphere's own frame (so you can see it spin),
// contacts as a yellow dot with a green normal. Lines further away are dimmer; nothing is hidden.
//
//   left-drag  orbit the camera          wheel  zoom
//   1 tower   2 pyramid   3 mixed rain   4 slope (friction and rolling)   5 tumble (zero gravity)
//   Space  shoot a ball from the camera   F  drop 10 more random bodies
//   W  warm starting on/off   - / =  fewer / more solver sweeps   G  gyroscopic term: implicit / explicit / off
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
constexpr size_t kMaxBodies = 400;  // the 3D world still tests every pair
constexpr Real kDt = 1.0f / 120.0f;

enum class Scene { Tower = 1, Pyramid, Rain, Slope, Tumble };
const char* scene_name(Scene s) {
    switch (s) {
        case Scene::Tower: return "tower";
        case Scene::Pyramid: return "pyramid";
        case Scene::Rain: return "mixed rain";
        case Scene::Slope: return "slope";
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
        const SDL_Color c = b.type == BodyType::Static ? SDL_Color{150, 154, 170, 255} : SDL_Color{120, 200, 255, 255};
        if (b.shape.type == Shape::Type::Box) paint.box(b, c);
        else paint.sphere(b, c);
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
    Scene scene = static_cast<Scene>(std::clamp(arg_scene, 1, 5));
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
                    case SDLK_1: case SDLK_2: case SDLK_3: case SDLK_4: case SDLK_5:
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
            SDL_RenderDebugText(renderer, 16, 12, "drag orbit | wheel zoom | 1-5 scene | Space shoot | F drop 10 | R reload | Tab pause | H hide text");
            SDL_RenderDebugTextFormat(renderer, 16, 28, "scene: %s | bodies %d | contact pairs %d | step %.2f ms%s", scene_name(scene),
                                      static_cast<int>(world.bodies.size()) - 1, static_cast<int>(world.contacts().size()), step_ms,
                                      paused ? " | PAUSED" : "");
            SDL_RenderDebugTextFormat(renderer, 16, 44, "W warm starting: %s | -/= sweeps: %d | G gyroscopic: %s | C contacts",
                                      world.solver.warm_starting ? "on" : "OFF", world.solver.iterations, gyro_name(world.gyroscopic));
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
