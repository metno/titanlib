cmake_minimum_required(VERSION 2.8)

# Change policy to allow a version to be set
cmake_policy(SET CMP0048 NEW)
project(titanlib VERSION 0.1.0)

# Add swig subdirectory
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wextra -Wno-unused-parameter -fPIC")
set(CMAKE_BUILD_TYPE Release)

# Add libtitanlib
include_directories(./include)
file(GLOB SOURCES "src/*.cpp")
add_library(titanlib SHARED ${SOURCES})
set_target_properties(titanlib PROPERTIES PUBLIC_HEADER "include/titanlib.h")
install(TARGETS titanlib
      EXPORT CMakeSwigTargets
      PUBLIC_HEADER DESTINATION include
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib
      )

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
add_subdirectory(SWIG)
