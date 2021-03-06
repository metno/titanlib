cmake_minimum_required(VERSION 2.8)

# Determine titanlib version
file (READ "${CMAKE_CURRENT_SOURCE_DIR}/include/titanlib.h" VERSION0)
project(titanlib)
string(REGEX REPLACE "^.+define TITANLIB_VERSION \"([^\n]+)\".+" "\\1" PROJECT_VERSION ${VERSION0})

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set compiler flags
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "CMAKE_BUILD_TYPE is not set, switching to RELEASE")
   set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
      "Choose the type of build, options are: NONE DEBUG RELEASE."
      FORCE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wfatal-errors -Wno-reorder -Wno-sign-compare -Wl,-init,initialize_omp")
if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
   set(CMAKE_CXX_FLAGS_DEBUG "-O0 -fPIC -g -pg -rdynamic --coverage")
   set(CMAKE_LFLAGS_DEBUG " --coverage")
   add_definitions("-DDEBUG")
elseif(CMAKE_BUILD_TYPE MATCHES "RELEASE")
   find_package(OpenMP)
   set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fPIC ${OpenMP_CXX_FLAGS}")
endif()

# Add libtitanlib
include_directories(./include)
file(GLOB SOURCES "src/*.cpp" "src/*.c")
add_library(titanlib ${SOURCES})
set_target_properties(titanlib PROPERTIES PUBLIC_HEADER "include/titanlib.h")
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
target_link_libraries(titanlib "${Boost_LIBRARIES}")
target_link_libraries(titanlib "${GSL_LIBRARIES}")
target_include_directories(titanlib PUBLIC "${Boost_INCLUDE_DIRS}")
target_include_directories(titanlib PUBLIC "${GSL_INCLUDE_DIRS}")
install(TARGETS titanlib
    EXPORT CMakeSwigTargets
    PUBLIC_HEADER DESTINATION include
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    )

# SWIG subdirectory
file(GLOB SWIG_INTERFACE "extras/SWIG/titanlib.i")
message('${SWIG_INTERFACE}')
add_subdirectory(extras/SWIG)
add_subdirectory(extras/python)

# Documentation
option(BUILD_DOC "Build documentation" ON)
find_package(Doxygen)
if (DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/config.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/docs/DoxygenLayout.xml" "DoxygenLayout.xml" @ONLY)
    add_custom_target(docs ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API
        documentation with Doxygen"
        VERBATIM )
else()
    message(STATUS "Documentation generation disabled. Doxygen needs to be installed to generate documentation.")
endif()

