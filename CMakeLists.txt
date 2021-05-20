cmake_minimum_required(VERSION 3.19.0)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}-std=c++1y")
project(DiplomProject)

#find_package(BLAS REQUIRED)
find_package(CUDA)
if (CUDA_FOUND)
    message("CUDA found")
else()
    message("CUDA not found, doing something alternatively")
endif()

file(GLOB headers "*.h*")
file(GLOB sources "*.cpp")

include_directories(#${BLAS_INCLUDE_DIRS}
                    ${CUDA_INCLUDE_DIRS}
                    ${CUDA_CUBLAS_DIRS})

add_executable(${PROJECT_NAME} ${headers} ${sources})

target_link_libraries(${PROJECT_NAME} PUBLIC ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusolver_LIBRARY}) #${BLAS_LIBRARIES})
#cuda_add_cublas_to_target(${PROJECT_NAME})