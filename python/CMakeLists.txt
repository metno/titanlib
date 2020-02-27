find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

# Files to install with Python
file(GLOB CPP_SOURCE_FILES ${CMAKE_SOURCE_DIR}/src)
file(GLOB CPP_INCLUDE_FILES ${CMAKE_SOURCE_DIR}/include)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CPP_SOURCE_FILES}
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CPP_INCLUDE_FILES}
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
 file(COPY ${SWIG_INTERFACE}
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# Configure setup.py and copy to output directory
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.in.py)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})
     #file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
     #        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Install target to call setup.py
add_custom_target(install-python
        COMMAND python3 ${SETUP_PY_OUT} install)

add_custom_target(install-python-user
        COMMAND python3 ${SETUP_PY_OUT} install --user)

add_custom_target(build-python
        COMMAND python3 ${SETUP_PY_OUT} build_ext)

add_custom_target(develop-python
        COMMAND python3 ${SETUP_PY_OUT} develop)

