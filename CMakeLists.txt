cmake_minimum_required(VERSION 3.25)
project(physics_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(-Wall -Wextra -Wconversion -Wshadow)

option(PHYS_SANITIZE "Build with ASan + UBSan" OFF)
if(PHYS_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

# Engine: pure library, no graphics dependency.
add_library(phys STATIC engine/src/shapes.cpp engine/src/body.cpp engine/src/collision.cpp engine/src/solver.cpp engine/src/dynamic_tree.cpp engine/src/broadphase.cpp engine/src/world.cpp)
target_include_directories(phys PUBLIC engine/include)

enable_testing()
add_subdirectory(tests)
add_subdirectory(bench)

find_package(SDL3 CONFIG QUIET)
if(SDL3_FOUND)
  add_subdirectory(demo)
else()
  message(STATUS "SDL3 not found, skipping demo")
endif()
