cmake_minimum_required(VERSION 4.0)
project(MES)

set(CMAKE_CXX_STANDARD 20)

add_executable(MES main.cpp
        source/data.cpp
        source/gaussLegendre.cpp)
