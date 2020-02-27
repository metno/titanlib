# Set up swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

option(BUILD_R "Build R SWIG module" ON)
if(BUILD_R)
    add_subdirectory(R)
endif()
