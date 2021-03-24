cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}-std=c++1y")
project(DiplomProject)

#set (BOOST_COMPONENTS numeric)
#set (Boost_USE_STATIC_LIBS ON)
#set (Boost_USE_MULTITHREADED ON)

file(GLOB headers "*.h*")
file(GLOB sources "*.cpp")

#find_package(Boost 1.74.0 REQUIRED) #COMPONENTS ${BOOST_COMPONENTS} REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${headers} ${sources})
#target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})