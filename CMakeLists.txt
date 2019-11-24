cmake_minimum_required(VERSION 3.10)
project(physics-sim)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS -pthread)

add_subdirectory(src)
