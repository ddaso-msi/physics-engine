add_executable(phys_tests main.cpp test_math.cpp test_integrate.cpp test_body.cpp test_collision.cpp test_world.cpp test_broadphase.cpp test_joints.cpp test_stage9.cpp)
target_link_libraries(phys_tests PRIVATE phys_scenes)
add_test(NAME phys_tests COMMAND phys_tests)
