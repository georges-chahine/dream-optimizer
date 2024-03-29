cmake_minimum_required(VERSION 2.8.3)
project(dream_offline_optimizer)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)
#set(CMAKE_CXX_STANDARD 11)
#add_compile_options(-std=c++11)
## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/thirdparty/g2o/cmake_modules ) 
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/thirdparty/g2o/EXTERNAL/csparse ) 
set (G2O_ROOT  ${PROJECT_SOURCE_DIR}/thirdparty/g2o/build)
set (CSPARSE_ROOT  ${PROJECT_SOURCE_DIR}/thirdparty/g2o/EXTERNAL/csparse)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rosbag
  tf
  cv_bridge
  tf_conversions
  sensor_msgs
  std_msgs
  pcl_ros

)
find_package(G2O REQUIRED)
FIND_PACKAGE(SuiteSparse)
#find_package(RTABMap)



#include_directories(include)
# flags
add_definitions("-DENABLE_SSE")
#set(CMAKE_CXX_FLAGS
#   "${SSE_FLAGS} -O3 -g -std=gnu++11 -march=native"
#   "${SSE_FLAGS} -O3 -g -std=c++0x -fno-omit-frame-pointer"
#)

#if (MSVC)
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
#endif (MSVC)


set(CMAKE_BUILD_TYPE Debug)

#add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/g2o)
link_directories(${PROJECT_SOURCE_DIR}/thirdparty/g2o/lib ${PROJECT_SOURCE_DIR}/thirdparty/g2o/EXTERNAL/csparse)
## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
 ${catkin_INCLUDE_DIRS} src ${PROJECT_SOURCE_DIR}/thirdparty/g2o ${PROJECT_SOURCE_DIR}/thirdparty/g2o/build ${PROJECT_SOURCE_DIR}/thirdparty/g2o/EXTERNAL/csparse ${SUITESPARSE_INCLUDE_DIRS} ${CSPARSE_INCLUDE_DIR}
)




## Declare a C++ library
# add_library(${PROJECT_NAME}
#   src/${PROJECT_NAME}/robot_tf.cpp
# )

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
## With catkin_make all packages are built within a single CMake context
## The recommended prefix ensures that target names across packages don't collide
# add_executable(${PROJECT_NAME}_node src/robot_tf_node.cpp)

SET (G2O_LIBS cxsparse g2o_csparse_extension g2o_cli g2o_core g2o_interface g2o_parser g2o_solver_cholmod
g2o_solver_csparse g2o_solver_dense g2o_solver_pcg g2o_stuff g2o_types_icp
g2o_types_sba g2o_types_sim3 g2o_types_slam2d g2o_types_slam3d)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/thirdparty/g2o/)
set(LIBRARY_OUTPUT_PATH)
add_executable(dream_offline_optimizer src/main.cpp)


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/thirdparty/g2o/bin)


####### g2o custom types ###############################
#ADD_SUBDIRECTORY(src/custom_types)

## Rename C++ executable without prefix
## The above recommended prefix causes long target names, the following renames the
## target back to the shorter version for ease of user use
## e.g. "rosrun someones_pkg node" instead of "rosrun someones_pkg someones_pkg_node"
# set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME node PREFIX "")

## Add cmake target dependencies of the executable
## same as for the library above
# add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
#message("libpointchecker_LIBRARIES: ${libpointchecker_LIBRARIES}")
## Specify libraries to link a library or executable target against
 target_link_libraries(dream_offline_optimizer yaml-cpp
   ${catkin_LIBRARIES} ${G2O_LIBS} cholmod #${RTABMap_LIBRARIES} #g2o_custom_types cholmod
 )

#############
## Install ##
#############

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# install(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables and/or libraries for installation
# install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_node
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )

add_definitions(-g -ggdb -O3 -Wall)
#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_robot_tf.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
