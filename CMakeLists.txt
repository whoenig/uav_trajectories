cmake_minimum_required(VERSION 3.1)

project(uavTrajectoryHelpers)

# Enable C++14 and warnings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost 1.58 REQUIRED COMPONENTS program_options filesystem)

# Additional include folders
include_directories(
  mav_trajectory_generation/mav_trajectory_generation/include
  ${EIGEN3_INCLUDE_DIR}
)

add_library(mav_trajectory_generation
  mav_trajectory_generation/mav_trajectory_generation/src/motion_defines.cpp
  mav_trajectory_generation/mav_trajectory_generation/src/polynomial.cpp
  mav_trajectory_generation/mav_trajectory_generation/src/rpoly.cpp
  mav_trajectory_generation/mav_trajectory_generation/src/segment.cpp
  mav_trajectory_generation/mav_trajectory_generation/src/timing.cpp
  mav_trajectory_generation/mav_trajectory_generation/src/trajectory.cpp
  mav_trajectory_generation/mav_trajectory_generation/src/vertex.cpp
)

add_executable(genTrajectory
  src/genTrajectory.cpp
)
target_link_libraries(genTrajectory
  ${Boost_LIBRARIES}
  mav_trajectory_generation
  glog
  nlopt
)
