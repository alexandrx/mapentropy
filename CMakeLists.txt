cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(MeanMapEntropy)

set(CMAKE_CXX_STANDARD 14) # C++14...
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PCL 1.7 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS} include)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable (mean_map_entropy src/mean_map_entropy.cpp src/map_entropy.cpp)
target_link_libraries (mean_map_entropy ${PCL_LIBRARIES})
