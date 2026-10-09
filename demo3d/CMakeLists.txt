# Wireframe demo for the 3D engine.
add_executable(demo3d main.cpp)
target_link_libraries(demo3d PRIVATE phys3d SDL3::SDL3)
