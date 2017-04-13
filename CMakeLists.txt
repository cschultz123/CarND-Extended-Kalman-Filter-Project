cmake_minimum_required (VERSION 3.5)

project(ExtendedKF)

add_definitions(-std=c++0x)

set(sources
    src/FusionEKF.cpp
    src/kalman_filter.cpp
    src/main.cpp
    src/tools.cpp)

# Trivial example using gtest and gmock
add_executable(ExtendedKF ${sources})