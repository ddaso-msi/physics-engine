#pragma once
// Solid, shaded drawing of the 3D world with nothing but SDL_RenderGeometry: every surface is cut into
// triangles, lit by one sun, and drawn back to front (the painter's algorithm), with faces turned away
// from the camera skipped. That is right for convex bodies that do not pass through each other, which is
// all the engine makes; large flat faces are cut into tiles so they sort sensibly against small bodies.
// Guides (lines, dots, labels) are drawn after, on top of everything.
#include "camera.hpp"
#include "hud.hpp"

#include <level.hpp>

#include <vector>

namespace puzzle3d::app {

struct Shade {
    SDL_Color base;
    Uint8 alpha = 255;
};

class Painter {
public:
    Painter(SDL_Renderer* ren, const demo3d::Camera& cam);

    void sky();
    void ground(const Body& ground);          // the top of the ground slab, checkered; drawn at once
    void shadow(const Body& b);               // a soft dark disc on the ground under a body; drawn at once
    void body(const Body& b, Shade shade);    // queued
    void zone_floor(Vec3 lo, Vec3 hi, SDL_Color c);  // a tinted rectangle at lo.y, queued
    void flush();                             // draws everything queued, far to near

    void line(Vec3 a, Vec3 b, SDL_Color c, bool dashed = false, float width = 2.0f);
    void marker(Vec3 p, Real radius, SDL_Color c);  // a round dot, sized by perspective
    void label(Vec3 p, const char* str, SDL_Color c);
    bool project(Vec3 p, SDL_FPoint& out, Real* depth = nullptr) const;

private:
    struct Tri {
        SDL_Vertex v[3];
        float depth;
    };
    void tri(Vec3 a, Vec3 b, Vec3 c, Vec3 na, Vec3 nb, Vec3 nc, Shade s);
    void quad(Vec3 a, Vec3 b, Vec3 c, Vec3 d, Vec3 n, Shade s, int tiles_u = 1, int tiles_v = 1);
    void round_body(const Body& b, Shade s);  // spheres and capsules
    SDL_FColor lit(SDL_Color c, Uint8 alpha, Vec3 normal) const;

    SDL_Renderer* ren_;
    const demo3d::Camera& cam_;
    Vec3 eye_;
    std::vector<Tri> tris_;
};

}  // namespace puzzle3d::app
