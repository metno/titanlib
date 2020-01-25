# Include python
find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIR})

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(../${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
#include_directories(../..)

# Add swig module
swig_add_library(rtitanlib TYPE SHARED LANGUAGE R SOURCES ../${PROJECT_NAME}.i)
swig_link_libraries(rtitanlib ${PROJECT_NAME} ${R_LIBRARIES})

# Files to install with Python
#set(PYTHON_INSTALL_FILES
#        ${CMAKE_CURRENT_BINARY_DIR}/titanlib.py
#        ${CMAKE_CURRENT_BINARY_DIR}/_titanlib.so)


# Declare install target for python
#install(TARGETS swig_example
#        COMMAND "${PYTHON_EXECUTABLE} setup.py"
#        COMPONENT swig-python)

# Install target to call setup.py
