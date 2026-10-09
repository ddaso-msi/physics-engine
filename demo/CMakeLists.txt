# Stage 2: three integrators on one spring.
add_executable(demo_integrators integrators.cpp)
target_link_libraries(demo_integrators PRIVATE phys SDL3::SDL3)

# Stage 3: body state, impulses at a point.
add_executable(demo_bodies bodies.cpp)
target_link_libraries(demo_bodies PRIVATE phys SDL3::SDL3)

# Stage 4: narrow-phase collision debug view (no response).
add_executable(demo_collision collision.cpp)
target_link_libraries(demo_collision PRIVATE phys SDL3::SDL3)

# The main demo: a sandbox plus a gallery of scenes. Grows with each stage.
add_executable(demo main.cpp)
target_link_libraries(demo PRIVATE phys_scenes SDL3::SDL3)
