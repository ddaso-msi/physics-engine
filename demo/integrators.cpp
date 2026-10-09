// Stage 2 demo: the same spring (x'' = -k x) integrated three ways, side by side.
// Explicit Euler gains energy and flies off; semi-implicit Euler and Verlet stay bounded.
// The bar under each row is current energy / starting energy (clamped to the bar's width).
// R resets, Esc quits.
#include <SDL3/SDL.h>
#include <phys/integrate.hpp>
#include <phys/timestep.hpp>

#include <cmath>

using namespace phys;

namespace {

constexpr int kWidth = 960, kHeight = 640;
constexpr float kPixelsPerMeter = 90.0f;
constexpr float kStiffness = 36.0f;  // omega = 6 rad/s; large enough that Euler's error shows in seconds

struct Row {
    const char* label;
    Integrator kind;
    SDL_Color color;
    Particle p;
};

double energy(const Particle& p) {
    return 0.5 * double(p.vel.length_sq()) + 0.5 * double(kStiffness) * double(p.pos.length_sq());
}

void reset(Row (&rows)[3]) {
    for (Row& r : rows) r.p = Particle{{1.0f, 0.0f}, {0.0f, 0.0f}};
}

void draw_circle(SDL_Renderer* ren, float cx, float cy, float radius) {
    constexpr int kSegments = 28;
    SDL_FPoint pts[kSegments + 1];
    for (int i = 0; i <= kSegments; ++i) {
        float a = 2.0f * kPi * static_cast<float>(i) / kSegments;
        pts[i] = {cx + radius * std::cos(a), cy + radius * std::sin(a)};
    }
    SDL_RenderLines(ren, pts, kSegments + 1);
}

}  // namespace

int main() {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    if (!SDL_CreateWindowAndRenderer("physics-engine: integrators", kWidth, kHeight, 0, &window,
                                     &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    Row rows[3] = {
        {"explicit Euler", Integrator::ExplicitEuler, {255, 110, 100, 255}, {}},
        {"semi-implicit Euler", Integrator::SemiImplicitEuler, {255, 210, 90, 255}, {}},
        {"velocity Verlet", Integrator::VelocityVerlet, {110, 220, 140, 255}, {}},
    };
    reset(rows);
    const double e0 = energy(rows[0].p);

    FixedTimestep timestep(1.0f / 60.0f);
    Uint64 last = SDL_GetTicksNS();
    bool running = true;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = false;
            if (e.type == SDL_EVENT_KEY_DOWN) {
                if (e.key.key == SDLK_ESCAPE) running = false;
                if (e.key.key == SDLK_R) reset(rows);
            }
        }

        Uint64 now = SDL_GetTicksNS();
        Real frame = static_cast<Real>(static_cast<double>(now - last) * 1e-9);
        last = now;
        timestep.advance(frame, [&](Real dt) {
            for (Row& r : rows)
                step(r.p, dt, r.kind, [](const Particle& p) { return Vec2{-kStiffness * p.pos.x, -kStiffness * p.pos.y}; });
        });

        SDL_SetRenderDrawColor(renderer, 20, 22, 28, 255);
        SDL_RenderClear(renderer);
        for (int i = 0; i < 3; ++i) {
            const Row& r = rows[i];
            const float y = 130.0f + 190.0f * static_cast<float>(i);
            const float cx = kWidth * 0.5f;
            SDL_SetRenderDrawColor(renderer, 70, 74, 90, 255);
            SDL_RenderLine(renderer, 0.0f, y, static_cast<float>(kWidth), y);          // rest line
            SDL_RenderLine(renderer, cx, y - 14.0f, cx, y + 14.0f);                    // equilibrium tick
            SDL_SetRenderDrawColor(renderer, r.color.r, r.color.g, r.color.b, 255);
            const float px = cx + r.p.pos.x * kPixelsPerMeter;
            SDL_RenderLine(renderer, cx, y, px, y);                                    // the "spring"
            draw_circle(renderer, px, y, 12.0f);

            const float bar = 360.0f * static_cast<float>(std::fmin(energy(r.p) / e0, 1.0) * 0.5);
            SDL_FRect box{40.0f, y + 40.0f, bar, 8.0f};
            SDL_RenderFillRect(renderer, &box);
            SDL_RenderDebugText(renderer, 40.0f, y - 60.0f, r.label);
        }
        SDL_SetRenderDrawColor(renderer, 150, 154, 170, 255);
        SDL_RenderDebugText(renderer, 40.0f, 16.0f, "same spring, same dt=1/60. Bars: energy / start (half = 1x, full = 2x+). R resets.");
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
