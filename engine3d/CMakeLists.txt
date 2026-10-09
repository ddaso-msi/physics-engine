# The 3D engine (namespace phys3d). Header-only for now. It borrows Real from the 2D engine's math.hpp,
# and later the dimension-free timestep and solver settings, so it links "phys" for those headers.
add_library(phys3d INTERFACE)
target_include_directories(phys3d INTERFACE include)
target_link_libraries(phys3d INTERFACE phys)
