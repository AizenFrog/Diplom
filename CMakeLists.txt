cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}-std=c++1y")
project(DiplomProject)

file(GLOB headers "*.h*")
file(GLOB sources "*.cpp")

add_executable(${PROJECT_NAME} ${headers} ${sources})