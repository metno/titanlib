cmake_minimum_required(VERSION 2.8)

# Change policy to allow a version to be set
cmake_policy(SET CMP0048 NEW)
project(titanlib VERSION 0.1.0)

# Add swig subdirectory
add_subdirectory(src)
