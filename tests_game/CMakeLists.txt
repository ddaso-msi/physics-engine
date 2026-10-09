# Reuses the engine tests' dependency-free runner (test.hpp + main.cpp).
add_executable(game_tests ${CMAKE_SOURCE_DIR}/tests/main.cpp test_objectives.cpp test_levels.cpp)
target_include_directories(game_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(game_tests PRIVATE phys_game)
add_test(NAME game_tests COMMAND game_tests)
