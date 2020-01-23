cmake_minimum_required(VERSION 2.8)

# Change policy to allow a version to be set
cmake_policy(SET CMP0048 NEW)
project(titanlib VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wextra -Wno-unused-parameter -fPIC")

## Add libtitanlib
add_library(titanlib titanlib.cpp)

# Set up public headers for liba
set_target_properties(titanlib PROPERTIES PUBLIC_HEADER "titanlib.h")

# Create main and link to liba
# add_library(${PROJECT_NAME} titanlib.cpp)
#target_link_libraries(${PROJECT_NAME} titanlib)

# Set up public headers for swig-example
# set_target_properties(titanlib PROPERTIES PUBLIC_HEADER "titanlib.h")

# Add install target
install(TARGETS ${PROJECT_NAME} titanlib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static
    PUBLIC_HEADER DESTINATION include)

# Add swig subdirectory
add_subdirectory(swig)
