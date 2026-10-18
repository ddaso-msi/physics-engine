# Game rules: levels, objectives, counters, ray picking. No SDL here, so the tests can play a level headless.
add_library(phys_game3d STATIC level.cpp raycast.cpp objectives.cpp pieces.cpp levels/levels.cpp levels/knock_down.cpp
            levels/bridge.cpp levels/pendulum.cpp levels/chain_reaction.cpp levels/tower.cpp)
target_include_directories(phys_game3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(phys_game3d PUBLIC phys3d)

# The game itself: window, input, drawing. The camera is the 3D demo's (tested in tests3d).
if(SDL3_FOUND)
  add_executable(puzzle_lab_3d app/main.cpp app/hud.cpp app/painter.cpp app/scene.cpp)
  target_include_directories(puzzle_lab_3d PRIVATE ${CMAKE_SOURCE_DIR}/demo3d)
  target_link_libraries(puzzle_lab_3d PRIVATE phys_game3d SDL3::SDL3)
endif()
