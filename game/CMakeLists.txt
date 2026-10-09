# Game rules: levels, objectives, counters. No SDL here, so the tests can play a level headless.
add_library(phys_game STATIC level.cpp objectives.cpp pieces.cpp levels/levels.cpp levels/knock_down.cpp levels/bridge.cpp
            levels/pendulum.cpp levels/chain_reaction.cpp levels/tower.cpp)
target_include_directories(phys_game PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(phys_game PUBLIC phys_scenes)

# The game itself: window, input, drawing.
if(SDL3_FOUND)
  add_executable(puzzle_lab app/main.cpp app/render.cpp app/scene.cpp)
  target_link_libraries(puzzle_lab PRIVATE phys_game SDL3::SDL3)
endif()
