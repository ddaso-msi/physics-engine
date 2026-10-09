// Stage 3 demo: rigid body state. Bodies float in zero gravity with no collisions yet.
//
// Drag from a point ON a body and release: that applies an impulse at that point, toward where you
// released. Push through the centre and the body just translates; push off-centre and it also
// spins, by r x j / I. Green line = velocity, yellow dot = centre of mass.
//
//   G toggles gravity   R resets   Esc quits
// Bodies wrap around the window edges (demo only; the engine has no walls yet).
#include <SDL3/SDL.h>
#include <phys/body.hpp>
#include <phys/timestep.hpp>

#include <cmath>
#include <vector>

using namespace phys;

namespace {

constexpr int kWidth = 960, kHeight = 640;
constexpr float kPixelsPerMeter = 60.0f;
constexpr float kWorldW = static_cast<float>(kWidth) / kPixelsPerMeter;
constexpr float kWorldH = static_cast<float>(kHeight) / kPixelsPerMeter;

SDL_FPoint to_screen(Vec2 w) { return {w.x * kPixelsPerMeter, (kWorldH - w.y) * kPixelsPerMeter}; }
Vec2 to_world(float sx, float sy) { return {sx / kPixelsPerMeter, kWorldH - sy / kPixelsPerMeter}; }

std::vector<Body> make_scene() {
    std::vector<Body> bodies;
    bodies.emplace_back(Shape::make_polygon(Polygon::box(1.0f, 0.5f)), Vec2{3.0f, 7.0f}, 0.3f);

    bodies.emplace_back(Shape::make_circle(0.8f), Vec2{8.0f, 5.5f}, 0.0f);

    const Vec2 tri[3] = {{0, 0}, {2.2f, 0}, {0.5f, 1.8f}};
    bodies.emplace_back(Shape::make_polygon(*Polygon::from_points(tri)), Vec2{12.5f, 7.0f}, 0.0f);

    // A long thin bar is hard to spin about its centre: large inertia relative to its mass.
    bodies.emplace_back(Shape::make_polygon(Polygon::box(2.2f, 0.18f)), Vec2{8.0f, 2.2f}, -0.2f);
    return bodies;
}

void draw_body(SDL_Renderer* ren, const Body& b) {
    SDL_SetRenderDrawColor(ren, 120, 200, 255, 255);
    if (b.shape.type == Shape::Type::Circle) {
        constexpr int kSegments = 36;
        SDL_FPoint pts[kSegments + 1];
        for (int i = 0; i <= kSegments; ++i) {
            float a = 2.0f * kPi * static_cast<float>(i) / kSegments;
            pts[i] = to_screen(b.pos + rotate(b.q, Vec2{b.shape.circle.radius * std::cos(a),
                                                        b.shape.circle.radius * std::sin(a)}));
        }
        SDL_RenderLines(ren, pts, kSegments + 1);
        // One spoke so rotation is visible on a symmetric shape.
        SDL_FPoint c = to_screen(b.pos), tip = to_screen(b.pos + b.q.x_axis() * b.shape.circle.radius);
        SDL_RenderLine(ren, c.x, c.y, tip.x, tip.y);
    } else {
        const Polygon& p = b.shape.polygon;
        SDL_FPoint pts[Polygon::kMaxVertices + 1];
        for (int i = 0; i <= p.count; ++i) pts[i] = to_screen(apply(b.transform(), p.vertices[i % p.count]));
        SDL_RenderLines(ren, pts, p.count + 1);
    }

    SDL_FPoint c = to_screen(b.pos);
    SDL_SetRenderDrawColor(ren, 110, 220, 140, 255);
    SDL_FPoint v = to_screen(b.pos + b.vel * 0.4f);
    SDL_RenderLine(ren, c.x, c.y, v.x, v.y);
    SDL_SetRenderDrawColor(ren, 255, 210, 90, 255);
    SDL_FRect dot{c.x - 3.0f, c.y - 3.0f, 6.0f, 6.0f};
    SDL_RenderFillRect(ren, &dot);

    SDL_SetRenderDrawColor(ren, 150, 154, 170, 255);
    SDL_RenderDebugTextFormat(ren, c.x + 12.0f, c.y + 12.0f, "m=%.2f I=%.2f w=%.2f", static_cast<double>(b.mass),
                              static_cast<double>(b.inertia), static_cast<double>(b.w));
}

void wrap(Body& b) {
    if (b.pos.x < 0) b.pos.x += kWorldW;
    if (b.pos.x > kWorldW) b.pos.x -= kWorldW;
    if (b.pos.y < 0) b.pos.y += kWorldH;
    if (b.pos.y > kWorldH) b.pos.y -= kWorldH;
}

}  // namespace

int main() {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    if (!SDL_CreateWindowAndRenderer("physics-engine: rigid bodies", kWidth, kHeight, 0, &window, &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    std::vector<Body> bodies = make_scene();
    bool gravity_on = false;

    // Drag state: which body, and where on it (in its own frame so the grab point rides along).
    int grabbed = -1;
    Vec2 grab_local;
    Vec2 mouse;

    FixedTimestep timestep(1.0f / 60.0f);
    Uint64 last = SDL_GetTicksNS();
    bool running = true;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = false;
            if (e.type == SDL_EVENT_KEY_DOWN) {
                if (e.key.key == SDLK_ESCAPE) running = false;
                if (e.key.key == SDLK_G) gravity_on = !gravity_on;
                if (e.key.key == SDLK_R) { bodies = make_scene(); grabbed = -1; }
            }
            if (e.type == SDL_EVENT_MOUSE_MOTION) mouse = to_world(e.motion.x, e.motion.y);
            if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN && e.button.button == SDL_BUTTON_LEFT) {
                mouse = to_world(e.button.x, e.button.y);
                for (int i = static_cast<int>(bodies.size()) - 1; i >= 0; --i) {
                    if (bodies[static_cast<size_t>(i)].contains(mouse)) {
                        grabbed = i;
                        grab_local = apply_inv(bodies[static_cast<size_t>(i)].transform(), mouse);
                        break;
                    }
                }
            }
            if (e.type == SDL_EVENT_MOUSE_BUTTON_UP && e.button.button == SDL_BUTTON_LEFT && grabbed >= 0) {
                mouse = to_world(e.button.x, e.button.y);
                Body& b = bodies[static_cast<size_t>(grabbed)];
                Vec2 hit = apply(b.transform(), grab_local);
                // Scale by mass so every body gets a comparable velocity kick for the same drag.
                b.apply_impulse_at((mouse - hit) * (b.mass * 2.0f), hit);
                grabbed = -1;
            }
        }

        Uint64 now = SDL_GetTicksNS();
        Real frame = static_cast<Real>(static_cast<double>(now - last) * 1e-9);
        last = now;
        const Vec2 g = gravity_on ? Vec2{0.0f, -9.81f} : Vec2{};
        timestep.advance(frame, [&](Real dt) {
            for (Body& b : bodies) {
                b.integrate(dt, g);
                wrap(b);
            }
        });

        SDL_SetRenderDrawColor(renderer, 20, 22, 28, 255);
        SDL_RenderClear(renderer);
        for (const Body& b : bodies) draw_body(renderer, b);

        if (grabbed >= 0) {
            const Body& b = bodies[static_cast<size_t>(grabbed)];
            SDL_FPoint a = to_screen(apply(b.transform(), grab_local)), m = to_screen(mouse);
            SDL_SetRenderDrawColor(renderer, 255, 110, 100, 255);
            SDL_RenderLine(renderer, a.x, a.y, m.x, m.y);
        }
        SDL_SetRenderDrawColor(renderer, 150, 154, 170, 255);
        SDL_RenderDebugText(renderer, 16.0f, 12.0f,
                            "drag from a body and release to push it at that point | G gravity | R reset");
        SDL_RenderDebugText(renderer, 16.0f, 28.0f, gravity_on ? "gravity: on" : "gravity: off");
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
