# Stage 2: three integrators on one spring.
add_executable(demo_integrators integrators.cpp)
target_link_libraries(demo_integrators PRIVATE phys SDL3::SDL3)

# Stage 3: body state, impulses at a point.
add_executable(demo_bodies bodies.cpp)
target_link_libraries(demo_bodies PRIVATE phys SDL3::SDL3)

# Latest stage (currently 4: collision detection debug view). Grows with each stage.
add_executable(demo main.cpp)
target_link_libraries(demo PRIVATE phys SDL3::SDL3)
