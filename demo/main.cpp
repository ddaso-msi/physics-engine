// Stage 4 demo: collision detection debug view. No response yet, so bodies pass through each other.
//
// Drag a body with the mouse; rotate the body under the cursor with the wheel or Q / E.
// Overlapping pairs turn red and show their manifold:
//   yellow squares = contact points (1 or 2)      green line = contact normal, pointing A -> B
//   text           = minimum translation depth and contact count
// Try: lay a box flat on another (2 points), stand it on a corner (1 point), roll a circle into a
// corner (vertex region), push one box through another to see the reference face switch.
//
//   R resets   Esc quits
#include <SDL3/SDL.h>
#include <phys/collision.hpp>

#include <cmath>
#include <vector>

using namespace phys;

namespace {

constexpr int kWidth = 960, kHeight = 640;
constexpr float kPixelsPerMeter = 60.0f;
constexpr float kWorldH = static_cast<float>(kHeight) / kPixelsPerMeter;

SDL_FPoint to_screen(Vec2 w) { return {w.x * kPixelsPerMeter, (kWorldH - w.y) * kPixelsPerMeter}; }
Vec2 to_world(float sx, float sy) { return {sx / kPixelsPerMeter, kWorldH - sy / kPixelsPerMeter}; }

std::vector<Body> make_scene() {
    std::vector<Body> bodies;
    bodies.emplace_back(Shape::make_polygon(Polygon::box(1.5f, 0.6f)), Vec2{3.5f, 7.0f}, 0.0f);
    bodies.emplace_back(Shape::make_polygon(Polygon::box(0.9f, 0.9f)), Vec2{7.5f, 7.0f}, 0.5f);

    std::vector<Vec2> hex;
    for (int k = 0; k < 6; ++k) {
        float a = 2.0f * kPi * static_cast<float>(k) / 6.0f;
        hex.push_back({1.1f * std::cos(a), 1.1f * std::sin(a)});
    }
    bodies.emplace_back(Shape::make_polygon(*Polygon::from_points(hex)), Vec2{12.0f, 7.0f}, 0.0f);

    const Vec2 tri[3] = {{0, 0}, {2.4f, 0}, {0.6f, 1.9f}};
    bodies.emplace_back(Shape::make_polygon(*Polygon::from_points(tri)), Vec2{4.0f, 3.5f}, 0.0f);
    bodies.emplace_back(Shape::make_circle(0.9f), Vec2{8.5f, 3.5f}, 0.0f);
    bodies.emplace_back(Shape::make_circle(0.5f), Vec2{11.0f, 3.5f}, 0.0f);
    return bodies;
}

void draw_body(SDL_Renderer* ren, const Body& b, SDL_Color color) {
    SDL_SetRenderDrawColor(ren, color.r, color.g, color.b, 255);
    if (b.shape.type == Shape::Type::Circle) {
        constexpr int kSegments = 36;
        SDL_FPoint pts[kSegments + 1];
        const float r = b.shape.circle.radius;
        for (int i = 0; i <= kSegments; ++i) {
            float a = 2.0f * kPi * static_cast<float>(i) / kSegments;
            pts[i] = to_screen(b.pos + Vec2{r * std::cos(a), r * std::sin(a)});
        }
        SDL_RenderLines(ren, pts, kSegments + 1);
        SDL_FPoint c = to_screen(b.pos), tip = to_screen(b.pos + b.q.x_axis() * r);
        SDL_RenderLine(ren, c.x, c.y, tip.x, tip.y);
    } else {
        const Polygon& p = b.shape.polygon;
        SDL_FPoint pts[Polygon::kMaxVertices + 1];
        for (int i = 0; i <= p.count; ++i) pts[i] = to_screen(apply(b.transform(), p.vertices[i % p.count]));
        SDL_RenderLines(ren, pts, p.count + 1);
    }
}

// Topmost body containing the point, or -1.
int pick(const std::vector<Body>& bodies, Vec2 p) {
    for (int i = static_cast<int>(bodies.size()) - 1; i >= 0; --i)
        if (bodies[static_cast<size_t>(i)].contains(p)) return i;
    return -1;
}

}  // namespace

int main() {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    if (!SDL_CreateWindowAndRenderer("physics-engine: collision detection", kWidth, kHeight, 0, &window,
                                     &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    std::vector<Body> bodies = make_scene();
    int grabbed = -1;
    Vec2 grab_offset;  // body position minus mouse position at the moment of the grab
    Vec2 mouse;

    bool running = true;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = false;
            if (e.type == SDL_EVENT_KEY_DOWN) {
                if (e.key.key == SDLK_ESCAPE) running = false;
                if (e.key.key == SDLK_R) { bodies = make_scene(); grabbed = -1; }
                if (e.key.key == SDLK_Q || e.key.key == SDLK_E) {
                    int target = grabbed >= 0 ? grabbed : pick(bodies, mouse);
                    if (target >= 0) {
                        Body& b = bodies[static_cast<size_t>(target)];
                        b.set_angle(b.angle + (e.key.key == SDLK_Q ? 0.06f : -0.06f));
                    }
                }
            }
            if (e.type == SDL_EVENT_MOUSE_MOTION) {
                mouse = to_world(e.motion.x, e.motion.y);
                if (grabbed >= 0) bodies[static_cast<size_t>(grabbed)].pos = mouse + grab_offset;
            }
            if (e.type == SDL_EVENT_MOUSE_WHEEL) {
                int target = grabbed >= 0 ? grabbed : pick(bodies, mouse);
                if (target >= 0) {
                    Body& b = bodies[static_cast<size_t>(target)];
                    b.set_angle(b.angle + e.wheel.y * 0.08f);
                }
            }
            if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN && e.button.button == SDL_BUTTON_LEFT) {
                mouse = to_world(e.button.x, e.button.y);
                grabbed = pick(bodies, mouse);
                if (grabbed >= 0) grab_offset = bodies[static_cast<size_t>(grabbed)].pos - mouse;
            }
            if (e.type == SDL_EVENT_MOUSE_BUTTON_UP && e.button.button == SDL_BUTTON_LEFT) grabbed = -1;
        }

        // Narrow phase on every pair. (Stage 7 replaces this O(n^2) loop with a broad phase.)
        struct Hit { Manifold m; };
        std::vector<Hit> hits;
        std::vector<bool> touching(bodies.size(), false);
        for (size_t i = 0; i < bodies.size(); ++i)
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                Manifold m;
                if (collide(bodies[i], bodies[j], m)) {
                    hits.push_back({m});
                    touching[i] = touching[j] = true;
                }
            }

        SDL_SetRenderDrawColor(renderer, 20, 22, 28, 255);
        SDL_RenderClear(renderer);
        for (size_t i = 0; i < bodies.size(); ++i)
            draw_body(renderer, bodies[i], touching[i] ? SDL_Color{255, 110, 100, 255} : SDL_Color{120, 200, 255, 255});

        for (const Hit& h : hits) {
            for (int k = 0; k < h.m.count; ++k) {
                const SDL_FPoint p = to_screen(h.m.points[k].point);
                const SDL_FPoint tip = to_screen(h.m.points[k].point + h.m.normal * 0.7f);
                SDL_SetRenderDrawColor(renderer, 110, 220, 140, 255);
                SDL_RenderLine(renderer, p.x, p.y, tip.x, tip.y);
                SDL_SetRenderDrawColor(renderer, 255, 210, 90, 255);
                SDL_FRect dot{p.x - 4.0f, p.y - 4.0f, 8.0f, 8.0f};
                SDL_RenderFillRect(renderer, &dot);
            }
            const SDL_FPoint label = to_screen(h.m.points[0].point);
            SDL_SetRenderDrawColor(renderer, 200, 204, 220, 255);
            SDL_RenderDebugTextFormat(renderer, label.x + 10.0f, label.y + 10.0f, "depth %.3f  %d pt",
                                      static_cast<double>(h.m.depth), h.m.count);
        }

        SDL_SetRenderDrawColor(renderer, 150, 154, 170, 255);
        SDL_RenderDebugText(renderer, 16.0f, 12.0f, "drag bodies | wheel or Q/E rotates the body under the cursor | R reset");
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
