# The 3D engine (namespace phys3d). It borrows Real from the 2D engine's math.hpp, and later the
# dimension-free timestep and solver settings, so it links "phys" for those headers.
add_library(phys3d STATIC src/body.cpp src/collision.cpp)
target_include_directories(phys3d PUBLIC include)
target_link_libraries(phys3d PUBLIC phys)
