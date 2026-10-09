# Gallery scenes (cradle, rag doll, car). A separate library so tests and the demo share them.
add_library(phys_scenes STATIC scenes.cpp)
target_include_directories(phys_scenes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(phys_scenes PUBLIC phys)
