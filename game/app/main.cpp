// Physics Puzzle Lab: a puzzle game on top of the 2D engine.
//
// Level 1, Knock It Down: move the mouse to aim (the dotted arc is where the ball will go), click to
// fire. Knock every block off the platform with three balls.
//
//   R reset    Tab pause    H help    C show contacts    N next (after solving)    Esc back to the menu
//
// Checking a build without a display: arguments are run in order as a script, then the program exits.
//   --level N        open level N
//   --mouse X Y      put the cursor at a world position (metres)
//   --fire X Y       click at a world position
//   --steps N        advance the simulation N fixed steps
//   --key NAME       press a key: reset, pause, help, contacts, next, escape
//   --capture FILE   draw the current frame and save it as a BMP
#include "render.hpp"

#include <SDL3/SDL.h>
#include <levels/knock_down.hpp>
#include <phys/timestep.hpp>

#include <cstdlib>
#include <cstring>

using namespace puzzle;
using namespace puzzle::app;
using phys::Vec2;

namespace {

constexpr float kBarHeight = 88.0f;
constexpr int kLevelCount = 5;
constexpr int kPlayableLevels = 1;  // grows as levels are added

struct LevelCard {
    const char* name;
    const char* objective;
};
constexpr LevelCard kCards[kLevelCount] = {
    {"Knock It Down", "Knock every block off the platform."},
    {"Bridge Builder", "Build a bridge that carries a load across the gap."},
    {"Pendulum Smash", "Swing a weight into the target."},
    {"Chain Reaction", "Set off a chain of collisions to reach the target."},
    {"Impossible Tower", "Build a tall tower that stays standing."},
};

enum class Screen { Menu, Play };
enum class Action { None, Reset, Pause, Help, Menu, Next, Contacts };

struct App {
    Screen screen = Screen::Menu;
    KnockDown level;
    int attempt = 1;
    bool paused = false;
    bool show_help = false;
    bool show_contacts = false;
    bool running = true;
    float mouse_x = 0, mouse_y = 0;
    phys::FixedTimestep timestep{kTimeStep};
};

constexpr Button kButtons[] = {
    {{600.0f, 10.0f, 160.0f, 36.0f}, "R Reset"},
    {{768.0f, 10.0f, 160.0f, 36.0f}, "Tab Pause"},
    {{936.0f, 10.0f, 160.0f, 36.0f}, "H Help"},
    {{1104.0f, 10.0f, 160.0f, 36.0f}, "Esc Menu"},
};
constexpr Action kButtonActions[] = {Action::Reset, Action::Pause, Action::Help, Action::Menu};

SDL_FRect card_rect(int i) { return {140.0f, 250.0f + 78.0f * static_cast<float>(i), 1000.0f, 64.0f}; }
bool inside(const SDL_FRect& r, float x, float y) { return x >= r.x && x < r.x + r.w && y >= r.y && y < r.y + r.h; }

void start_level(App& app, int index) {
    if (index < 0 || index >= kPlayableLevels) return;
    app.level.reset();
    app.attempt = 1;
    app.paused = app.show_help = false;
    app.timestep = phys::FixedTimestep(kTimeStep);
    app.screen = Screen::Play;
}

void apply(App& app, Action action) {
    switch (action) {
        case Action::Reset:
            if (app.level.shots_used() > 0) ++app.attempt;  // resetting an untouched level is not a new attempt
            app.level.reset();
            app.paused = false;
            app.timestep = phys::FixedTimestep(kTimeStep);
            break;
        case Action::Pause: app.paused = !app.paused; break;
        case Action::Help: app.show_help = !app.show_help; break;
        case Action::Contacts: app.show_contacts = !app.show_contacts; break;
        case Action::Next:  // only once solved; with one level so far, "next" is the level list
            if (app.level.status == Status::Won) app.screen = Screen::Menu;
            break;
        case Action::Menu:
            if (app.show_help) app.show_help = false;  // Esc closes the overlay first
            else app.screen = Screen::Menu;
            break;
        case Action::None: break;
    }
}

void click(App& app, float x, float y) {
    if (app.screen == Screen::Menu) {
        for (int i = 0; i < kLevelCount; ++i)
            if (inside(card_rect(i), x, y)) start_level(app, i);
        return;
    }
    if (app.show_help) {  // any click dismisses the overlay and does nothing else
        app.show_help = false;
        return;
    }
    for (size_t i = 0; i < std::size(kButtons); ++i)
        if (kButtons[i].hit(x, y)) {
            apply(app, kButtonActions[i]);
            return;
        }
    if (y > kBarHeight && !app.paused) app.level.fire(to_world(x, y));
}

void key(App& app, SDL_Keycode k) {
    if (app.screen == Screen::Menu) {
        if (k == SDLK_Q) app.running = false;
        if (k == SDLK_RETURN) start_level(app, 0);
        if (k >= SDLK_1 && k <= SDLK_5) start_level(app, static_cast<int>(k - SDLK_1));
        return;
    }
    switch (k) {
        case SDLK_R: apply(app, Action::Reset); break;
        case SDLK_TAB: apply(app, Action::Pause); break;
        case SDLK_H: apply(app, Action::Help); break;
        case SDLK_C: apply(app, Action::Contacts); break;
        case SDLK_N: apply(app, Action::Next); break;
        case SDLK_ESCAPE: apply(app, Action::Menu); break;
        default: break;
    }
}

void handle_event(App& app, const SDL_Event& e) {
    switch (e.type) {
        case SDL_EVENT_QUIT: app.running = false; break;
        case SDL_EVENT_KEY_DOWN:
            if (!e.key.repeat) key(app, e.key.key);
            break;
        case SDL_EVENT_MOUSE_MOTION:
            app.mouse_x = e.motion.x;
            app.mouse_y = e.motion.y;
            break;
        case SDL_EVENT_MOUSE_BUTTON_DOWN:  // one event per press: holding the button does not repeat
            app.mouse_x = e.button.x;
            app.mouse_y = e.button.y;
            if (e.button.button == SDL_BUTTON_LEFT) click(app, e.button.x, e.button.y);
            break;
        default: break;
    }
}

void update(App& app, phys::Real frame_seconds) {
    // Paused, in the menu or reading the help: bank no time, so play resumes exactly where it stopped.
    if (app.screen != Screen::Play || app.paused || app.show_help) return;
    app.timestep.advance(frame_seconds, [&](phys::Real dt) { app.level.step(dt); });
}

// ---- drawing ----

void draw_grid(SDL_Renderer* ren) {
    set_color(ren, color::kGrid);
    for (int x = 1; x < static_cast<int>(kWorldW); ++x) {
        const SDL_FPoint p = to_screen({static_cast<float>(x), 0});
        SDL_RenderLine(ren, p.x, 0, p.x, static_cast<float>(kHeight));
    }
    for (int y = 1; y < static_cast<int>(kWorldH) + 1; ++y) {
        const SDL_FPoint p = to_screen({0, static_cast<float>(y)});
        SDL_RenderLine(ren, 0, p.y, static_cast<float>(kWidth), p.y);
    }
}

void draw_aim(SDL_Renderer* ren, const KnockDown& level, Vec2 aim) {
    // The same arithmetic the engine uses for a body in free flight, so the dots are where the ball goes.
    Vec2 p = KnockDown::kMuzzle, v = level.launch_velocity(aim);
    for (int i = 1; i <= 150; ++i) {
        v += level.world.gravity * kTimeStep;
        p += v * kTimeStep;
        if (p.y < kFloorTop || p.x < 0 || p.x > kWorldW) break;
        if (i % 5 != 0) continue;
        const Uint8 fade = static_cast<Uint8>(230 - i);
        fill_disc(ren, to_screen(p), 2.5f, {fade, fade, fade, 255});
    }
    const Vec2 dir = v.length_sq() > 0 ? level.launch_velocity(aim).normalized() : Vec2{1, 0};
    const SDL_FPoint a = to_screen(KnockDown::kMuzzle), b = to_screen(KnockDown::kMuzzle + dir * 0.75f);
    set_color(ren, color::kToolEdge);
    for (float o = -1.0f; o <= 1.0f; o += 1.0f) SDL_RenderLine(ren, a.x + o * dir.y, a.y + o * dir.x, b.x + o * dir.y, b.y + o * dir.x);
    text(ren, a.x - 44.0f, a.y - 58.0f, 2.0f, color::kToolEdge, "%2.0f m/s", static_cast<double>(level.launch_velocity(aim).length()));
}

void draw_world(SDL_Renderer* ren, const App& app) {
    const KnockDown& level = app.level;
    draw_grid(ren);

    const float px = KnockDown::kPlatformCenter.x;
    dashed_line(ren, px - 2.4f, px + 2.4f, KnockDown::kFallLine, color::kTarget);
    const SDL_FPoint label = to_screen({px + 2.5f, KnockDown::kFallLine});
    text(ren, label.x, label.y - 8.0f, 2.0f, color::kTarget, "fall line");

    size_t next_target = 0;
    for (size_t i = 0; i < level.world.bodies.size(); ++i) {
        const phys::Body& b = level.world.bodies[i];
        const bool is_target = next_target < level.targets.size() && static_cast<size_t>(level.targets[next_target]) == i;
        if (is_target) ++next_target;
        if (b.type == phys::BodyType::Static) draw_body(ren, b, color::kStatic, color::kStaticEdge);
        else if (is_target && b.pos.y < KnockDown::kFallLine) draw_body(ren, b, color::kDone, color::kDoneEdge);
        else if (is_target) draw_body(ren, b, color::kTarget, color::kTargetEdge);
        else draw_body(ren, b, color::kTool, color::kToolEdge);
    }

    // The launcher, with the next ball loaded while there is one.
    const SDL_FPoint muzzle = to_screen(KnockDown::kMuzzle);
    fill_disc(ren, muzzle, 0.36f * kPixelsPerMeter, color::kStaticEdge);
    fill_disc(ren, muzzle, 0.31f * kPixelsPerMeter, color::kStatic);
    if (level.shots_left > 0 && level.status == Status::Playing)
        fill_disc(ren, muzzle, KnockDown::kBallRadius * kPixelsPerMeter, level.can_fire() ? color::kTool : color::kDim);
    if (level.can_fire() && !app.paused && !app.show_help && app.mouse_y > kBarHeight)
        draw_aim(ren, level, to_world(app.mouse_x, app.mouse_y));

    if (app.show_contacts)
        for (const phys::ContactPair& c : level.world.contacts())
            for (int k = 0; k < c.manifold.count; ++k) {
                const SDL_FPoint p = to_screen(c.manifold.points[k].point);
                const SDL_FPoint tip = to_screen(c.manifold.points[k].point + c.manifold.normal * 0.3f);
                set_color(ren, color::kGood);
                SDL_RenderLine(ren, p.x, p.y, tip.x, tip.y);
                fill_rect(ren, {p.x - 2.0f, p.y - 2.0f, 4.0f, 4.0f}, {255, 240, 150, 255});
            }
}

void draw_panel(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color edge) {
    fill_rect(ren, r, color::kPanel);
    outline_rect(ren, r, edge);
    outline_rect(ren, {r.x + 1, r.y + 1, r.w - 2, r.h - 2}, edge);
}

void draw_hud(SDL_Renderer* ren, const App& app) {
    const KnockDown& level = app.level;
    fill_rect(ren, {0, 0, static_cast<float>(kWidth), kBarHeight}, color::kPanel);
    set_color(ren, color::kPanelEdge);
    SDL_RenderLine(ren, 0, kBarHeight, static_cast<float>(kWidth), kBarHeight);

    text(ren, 20.0f, 14.0f, 3.0f, color::kText, "1 %s", kCards[0].name);
    // What the game is doing right now, so aiming, paused and finished never look alike.
    const char* mode = app.paused ? "PAUSED" : level.status == Status::Won ? "SOLVED" : level.status == Status::Failed ? "FAILED" : level.can_fire() ? "AIM" : level.shots_left > 0 ? "RELOADING" : "WATCH";
    const SDL_Color mode_color = app.paused ? color::kTarget : level.status == Status::Won ? color::kGood : level.status == Status::Failed ? color::kBad : color::kToolEdge;
    text(ren, 400.0f, 18.0f, 2.0f, mode_color, "%s", mode);

    for (size_t i = 0; i < std::size(kButtons); ++i) draw_button(ren, kButtons[i], kButtons[i].hit(app.mouse_x, app.mouse_y));

    text(ren, 20.0f, 60.0f, 2.0f, color::kDim, "%s", kCards[0].objective);
    float x = 640.0f;
    text(ren, x, 60.0f, 2.0f, color::kDim, "Balls");
    for (int i = 0; i < KnockDown::kShots; ++i)
        fill_disc(ren, {x + 100.0f + 22.0f * static_cast<float>(i), 68.0f}, 8.0f, i < level.shots_left ? color::kTool : SDL_Color{50, 56, 74, 255});
    x += 190.0f;
    text(ren, x, 60.0f, 2.0f, level.fallen > 0 ? color::kTarget : color::kDim, "Down %d/%d", level.fallen, static_cast<int>(level.targets.size()));
    x += 150.0f;
    text(ren, x, 60.0f, 2.0f, color::kDim, "Try %d", app.attempt);
    x += 120.0f;
    text(ren, x, 60.0f, 2.0f, color::kDim, "%5.1f s", level.sim_time);

    const float cx = static_cast<float>(kWidth) * 0.5f;
    if (app.show_help) {
        // the overlay below covers the result banner
    } else if (level.status == Status::Won) {
        const SDL_FRect r{cx - 340.0f, 112.0f, 680.0f, 124.0f};
        draw_panel(ren, r, color::kGood);
        text_centered(ren, cx, r.y + 14.0f, 4.0f, color::kGood, "SOLVED!");
        char line[96];
        SDL_snprintf(line, sizeof(line), "%d ball%s, %.1f s, try %d", level.shots_used(), level.shots_used() == 1 ? "" : "s", level.win_time, app.attempt);
        text_centered(ren, cx, r.y + 60.0f, 2.0f, color::kText, line);
        text_centered(ren, cx, r.y + 92.0f, 2.0f, color::kDim, "N continue   R play again   Esc menu");
    } else if (level.status == Status::Failed) {
        const SDL_FRect r{cx - 300.0f, 112.0f, 600.0f, 112.0f};
        draw_panel(ren, r, color::kBad);
        text_centered(ren, cx, r.y + 14.0f, 3.0f, color::kBad, "OUT OF BALLS");
        char line[96];
        SDL_snprintf(line, sizeof(line), "%d of %d blocks down", level.fallen, static_cast<int>(level.targets.size()));
        text_centered(ren, cx, r.y + 50.0f, 2.0f, color::kText, line);
        text_centered(ren, cx, r.y + 80.0f, 2.0f, color::kDim, "R try again   Esc menu");
    }

    if (app.show_help) {
        const SDL_FRect r{cx - 420.0f, 150.0f, 840.0f, 400.0f};
        draw_panel(ren, r, color::kToolEdge);
        text_centered(ren, cx, r.y + 20.0f, 3.0f, color::kText, "HOW TO PLAY");
        const char* lines[] = {
            "Move the mouse to aim. The dots show the path.",
            "Aim further from the launcher for more speed.",
            "Click to fire. You have three balls.",
            "Every block must drop below the fall line.",
            "",
            "R      reset the level",
            "Tab    pause or resume",
            "N      continue after solving",
            "C      show contact points",
            "Esc    close this, or back to the menu",
        };
        float y = r.y + 72.0f;
        for (const char* l : lines) {
            text(ren, r.x + 40.0f, y, 2.0f, color::kText, "%s", l);
            y += 30.0f;
        }
    }
}

void draw_menu(SDL_Renderer* ren, const App& app) {
    draw_grid(ren);
    const float cx = static_cast<float>(kWidth) * 0.5f;
    text_centered(ren, cx, 70.0f, 6.0f, color::kText, "PHYSICS PUZZLE LAB");
    text_centered(ren, cx, 140.0f, 2.0f, color::kDim, "Five problems. Real physics. Experiment until it works.");
    text_centered(ren, cx, 200.0f, 2.0f, color::kToolEdge, "Pick a level");
    for (int i = 0; i < kLevelCount; ++i) {
        const SDL_FRect r = card_rect(i);
        const bool playable = i < kPlayableLevels;
        const bool hovered = playable && inside(r, app.mouse_x, app.mouse_y);
        fill_rect(ren, r, hovered ? SDL_Color{48, 56, 78, 255} : SDL_Color{33, 38, 52, 255});
        outline_rect(ren, r, hovered ? color::kToolEdge : color::kPanelEdge);
        text(ren, r.x + 20.0f, r.y + 20.0f, 3.0f, playable ? color::kTarget : color::kDim, "%d", i + 1);
        text(ren, r.x + 70.0f, r.y + 12.0f, 2.0f, playable ? color::kText : color::kDim, "%s", kCards[i].name);
        text(ren, r.x + 70.0f, r.y + 38.0f, 2.0f, color::kDim, "%s", playable ? kCards[i].objective : "Not built yet.");
        if (playable) text(ren, r.x + r.w - 90.0f, r.y + 24.0f, 2.0f, color::kGood, "PLAY");
    }
    text_centered(ren, cx, 660.0f, 2.0f, color::kDim, "Click a level or press its number.   Q quit");
}

void draw(SDL_Renderer* ren, const App& app) {
    set_color(ren, color::kBackground);
    SDL_RenderClear(ren);
    if (app.screen == Screen::Menu) {
        draw_menu(ren, app);
    } else {
        draw_world(ren, app);
        draw_hud(ren, app);
    }
}

// ---- scripted run (see the top of the file) ----

bool run_script(App& app, SDL_Renderer* ren, int argc, char** argv) {
    auto number = [&](int i) { return i < argc ? static_cast<float>(std::atof(argv[i])) : 0.0f; };
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!std::strcmp(a, "--level")) {
            start_level(app, static_cast<int>(number(++i)) - 1);
        } else if (!std::strcmp(a, "--mouse") || !std::strcmp(a, "--fire")) {
            const SDL_FPoint p = to_screen({number(i + 1), number(i + 2)});
            i += 2;
            app.mouse_x = p.x;
            app.mouse_y = p.y;
            if (!std::strcmp(a, "--fire")) click(app, p.x, p.y);
        } else if (!std::strcmp(a, "--steps")) {
            const int n = static_cast<int>(number(++i));
            for (int s = 0; s < n; ++s) update(app, kTimeStep * 1.0001f);
        } else if (!std::strcmp(a, "--key") && i + 1 < argc) {
            const char* k = argv[++i];
            key(app, !std::strcmp(k, "reset") ? SDLK_R : !std::strcmp(k, "pause") ? SDLK_TAB : !std::strcmp(k, "help") ? SDLK_H
                     : !std::strcmp(k, "contacts") ? SDLK_C : !std::strcmp(k, "next") ? SDLK_N : SDLK_ESCAPE);
        } else if (!std::strcmp(a, "--capture") && i + 1 < argc) {
            draw(ren, app);
            SDL_Surface* shot = SDL_RenderReadPixels(ren, nullptr);
            const bool ok = shot && SDL_SaveBMP(shot, argv[++i]);
            SDL_DestroySurface(shot);
            if (!ok) {
                SDL_Log("capture failed: %s", SDL_GetError());
                return false;
            }
        } else {
            SDL_Log("unknown argument: %s", a);
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    if (!SDL_CreateWindowAndRenderer("Physics Puzzle Lab", kWidth, kHeight, 0, &window, &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    App app;
    int result = 0;
    if (argc > 1) {
        result = run_script(app, renderer, argc, argv) ? 0 : 1;
    } else {
        Uint64 last = SDL_GetTicksNS();
        while (app.running) {
            SDL_Event e;
            while (SDL_PollEvent(&e)) handle_event(app, e);
            const Uint64 now = SDL_GetTicksNS();
            update(app, static_cast<phys::Real>(static_cast<double>(now - last) * 1e-9));
            last = now;
            draw(renderer, app);
            SDL_RenderPresent(renderer);
        }
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return result;
}
