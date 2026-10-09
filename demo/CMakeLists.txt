# Stage 2: three integrators on one spring.
add_executable(demo_integrators integrators.cpp)
target_link_libraries(demo_integrators PRIVATE phys SDL3::SDL3)

# Stage 3+: rigid bodies. Grows with each stage.
add_executable(demo main.cpp)
target_link_libraries(demo PRIVATE phys SDL3::SDL3)
