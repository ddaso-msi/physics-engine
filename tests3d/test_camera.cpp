#include "test.hpp"
#include <camera.hpp>

using demo3d::Camera;
using demo3d::ScreenPoint;
using phys3d::Vec3;

#define CHECK_VEC(a, b, eps)                         \
    do {                                             \
        const Vec3 va_ = (a), vb_ = (b);             \
        CHECK_NEAR(va_.x, vb_.x, eps);               \
        CHECK_NEAR(va_.y, vb_.y, eps);               \
        CHECK_NEAR(va_.z, vb_.z, eps);               \
    } while (0)

TEST(camera_frame_is_orthonormal_and_looks_at_the_target) {
    Camera cam;
    cam.yaw = 0;
    cam.pitch = 0;
    cam.target = {0, 0, 0};
    cam.distance = 10;
    CHECK_VEC(cam.eye(), Vec3(0, 0, 10), 1e-4);     // yaw 0, pitch 0: straight back along +z
    CHECK_VEC(cam.forward(), Vec3(0, 0, -1), 1e-4);
    CHECK_VEC(cam.right(), Vec3(1, 0, 0), 1e-4);    // world +x is to the right
    CHECK_VEC(cam.up(), Vec3(0, 1, 0), 1e-4);

    for (float yaw : {-2.0f, 0.3f, 1.9f})
        for (float pitch : {-1.0f, 0.0f, 0.8f}) {
            Camera c;
            c.yaw = yaw;
            c.pitch = pitch;
            c.target = {1, 2, 3};
            CHECK_NEAR(distance(c.eye(), c.target), c.distance, 1e-3);
            CHECK_NEAR(dot(c.right(), c.forward()), 0, 1e-4);
            CHECK_NEAR(dot(c.up(), c.forward()), 0, 1e-4);
            CHECK_NEAR(c.up().length(), 1, 1e-4);
            CHECK_VEC(cross(c.right(), c.up()), -c.forward(), 1e-4);  // right-handed, looking down -z of its own frame
            CHECK(c.up().y > 0);                                       // never upside down
            CHECK_NEAR(c.to_view(c.target).z, c.distance, 1e-3);       // the target is straight ahead
        }
}

TEST(projection_centres_the_target_and_shrinks_with_distance) {
    Camera cam;
    const ScreenPoint centre = cam.to_screen(cam.to_view(cam.target));
    CHECK_NEAR(centre.x, cam.width / 2.0, 1e-2);
    CHECK_NEAR(centre.y, cam.height / 2.0, 1e-2);

    // To the camera's right is +x on screen; above is SMALLER y, because screen y runs downward.
    const ScreenPoint right = cam.to_screen(cam.to_view(cam.target + cam.right()));
    const ScreenPoint above = cam.to_screen(cam.to_view(cam.target + cam.up()));
    CHECK(right.x > centre.x + 1 && std::fabs(right.y - centre.y) < 0.01f);
    CHECK(above.y < centre.y - 1 && std::fabs(above.x - centre.x) < 0.01f);

    // Perspective: the same sideways offset twice as far away covers half as many pixels.
    const float near_px = cam.to_screen({1, 0, 5}).x - cam.width / 2.0f, far_px = cam.to_screen({1, 0, 10}).x - cam.width / 2.0f;
    CHECK_NEAR(near_px, 2 * far_px, 1e-2);

    // The vertical field of view fills the window: a point at tan(fov/2) above the axis lands on the top edge.
    const float edge = std::tan(cam.fov_y * 0.5f);
    CHECK_NEAR(cam.to_screen({0, edge * 7, 7}).y, 0, 1e-2);
    CHECK_NEAR(cam.to_screen({0, -edge * 7, 7}).y, cam.height, 1e-2);
}

TEST(segments_are_clipped_at_the_near_plane) {
    Camera cam;
    cam.yaw = 0;
    cam.pitch = 0;
    cam.target = {0, 0, 0};
    cam.distance = 10;  // eye at z = 10 looking toward -z
    ScreenPoint a, b;
    phys3d::Real depth = 0;

    CHECK(cam.project_segment({-1, 0, 0}, {1, 0, 0}, a, b, &depth));  // wholly in front
    CHECK_NEAR(depth, 10, 1e-3);
    CHECK(a.x < b.x);

    CHECK(!cam.project_segment({0, 0, 11}, {1, 0, 20}, a, b));  // wholly behind the camera

    // One end in front, one far behind: the visible part ends on the near plane, with finite coordinates.
    CHECK(cam.project_segment({1, 0, 0}, {1, 0, 30}, a, b, &depth));
    CHECK(std::isfinite(a.x) && std::isfinite(b.x) && std::isfinite(b.y));
    CHECK_NEAR(depth, (10 + cam.near_plane) / 2, 1e-3);
    // The cut end is the same world line seen very close up: far out to the side, not flipped to the left.
    CHECK(b.x > a.x);
}

TEST(orbit_and_zoom_stay_within_bounds) {
    Camera cam;
    cam.orbit(0, 10);
    CHECK(cam.pitch <= 1.5f && cam.up().y > 0);
    cam.orbit(0, -20);
    CHECK(cam.pitch >= -1.5f);
    for (int i = 0; i < 100; ++i) cam.zoom(0.5f);
    CHECK(cam.distance >= 2);
    for (int i = 0; i < 100; ++i) cam.zoom(2.0f);
    CHECK(cam.distance <= 80);
}
