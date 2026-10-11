# Not a test: prints timings. Build the release preset for meaningful numbers.
add_executable(phys3d_bench main.cpp)
target_link_libraries(phys3d_bench PRIVATE phys3d)
