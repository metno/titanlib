set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wextra -Wno-unused-parameter -fPIC")

## Add libtitanlib
include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(titanlib ${SOURCES})

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
add_subdirectory(swig)
