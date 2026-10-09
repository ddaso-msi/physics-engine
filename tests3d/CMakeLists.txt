# Reuses the 2D tests' dependency-free runner (test.hpp + main.cpp).
add_executable(phys3d_tests ${CMAKE_SOURCE_DIR}/tests/main.cpp test_math3d.cpp test_inertia3d.cpp test_body3d.cpp test_collision3d.cpp test_world3d.cpp test_camera.cpp)
target_include_directories(phys3d_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests ${CMAKE_SOURCE_DIR}/demo3d)
target_link_libraries(phys3d_tests PRIVATE phys3d)
add_test(NAME phys3d_tests COMMAND phys3d_tests)
