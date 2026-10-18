# Reuses the engine tests' dependency-free runner (test.hpp + main.cpp).
add_executable(game3d_tests ${CMAKE_SOURCE_DIR}/tests/main.cpp test_basics.cpp test_levels.cpp)
target_include_directories(game3d_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(game3d_tests PRIVATE phys_game3d)
add_test(NAME game3d_tests COMMAND game3d_tests)
