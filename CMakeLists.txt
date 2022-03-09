cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Default build type: Release" FORCE)
endif()
project(mean_map_entropy)

set(CMAKE_CXX_STANDARD 14) # C++14...
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PCL 1.8 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS} include)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(catkin QUIET)
#if catkin found, build ros version
if(catkin_FOUND)
 find_package(catkin REQUIRED COMPONENTS
   sensor_msgs
   pcl_ros
   pcl_conversions
   roscpp
)

 include_directories(
   include
   ${PCL_INCLUDE_DIRS}
   ${catkin_INCLUDE_DIRS}
 )

 catkin_package(
   INCLUDE_DIRS include
   CATKIN_DEPENDS 
     roscpp 
     pcl_ros
     std_msgs 
     sensor_msgs
   DEPENDS 
      PCL
 )

 add_executable (mean_map_entropy src/mean_map_entropy_node.cpp src/map_entropy.cpp)
 target_link_libraries (mean_map_entropy ${catkin_LIBRARIES})

 install(TARGETS mean_map_entropy
             ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
             LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
             RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
             )

 install(DIRECTORY launch/
         DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/launch
         PATTERN ".svn" EXCLUDE
         )
endif()

add_executable (mean_map_entropy_pcl src/mean_map_entropy.cpp src/map_entropy.cpp)
target_link_libraries (mean_map_entropy_pcl ${PCL_LIBRARIES})

 install(TARGETS mean_map_entropy_pcl
             ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
             LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
             RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
             )