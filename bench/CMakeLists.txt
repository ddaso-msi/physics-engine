# Not a test: prints timings. Build the release preset for meaningful numbers.
add_executable(phys_bench main.cpp)
target_link_libraries(phys_bench PRIVATE phys)
