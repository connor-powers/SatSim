cmake_minimum_required(VERSION 3.10)
project(SatSim VERSION 1.0.0)
set(CMAKE_DEBUG_POSTFIX d)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "--coverage -g -O0")
IF(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
ENDIF(NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
include_directories(include)

include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
)


FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4
)

FetchContent_MakeAvailable(json googletest Eigen)

add_library(SatSim SHARED src/utils.cpp src/Satellite.cpp src/PhasedArrayGroundStation.cpp)
add_library(SatSim::SatSim ALIAS SatSim)
target_link_libraries(SatSim PRIVATE nlohmann_json::nlohmann_json Eigen3::Eigen)
install(TARGETS SatSim DESTINATION lib)
install(FILES include/utils.h include/Satellite.h include/PhasedArrayGroundStation.h DESTINATION include)

add_executable(Satellite_tests tests/Satellite_tests.cpp)
add_executable(utils_tests tests/utils_tests.cpp)
add_executable(gs_tests tests/PhasedArrayGroundStation_tests.cpp)

target_link_libraries(Satellite_tests PRIVATE gtest_main SatSim::SatSim nlohmann_json::nlohmann_json Eigen3::Eigen)
target_link_libraries(utils_tests PRIVATE gtest_main SatSim::SatSim nlohmann_json::nlohmann_json Eigen3::Eigen)
target_link_libraries(gs_tests PRIVATE gtest_main SatSim::SatSim nlohmann_json::nlohmann_json Eigen3::Eigen)

