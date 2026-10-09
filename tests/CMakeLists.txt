add_executable(phys_tests main.cpp test_math.cpp test_integrate.cpp test_body.cpp test_collision.cpp test_world.cpp)
target_link_libraries(phys_tests PRIVATE phys)
add_test(NAME phys_tests COMMAND phys_tests)
