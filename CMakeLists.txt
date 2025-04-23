cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "--coverage -g -O0")
IF(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
ENDIF(NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
include_directories(include)

project(orbital_sim)
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



add_executable(run simulation_setup.cpp src/utils.cpp src/Satellite.cpp src/PhasedArrayGroundStation.cpp)
add_executable(Satellite_tests tests/Satellite_tests.cpp src/Satellite.cpp src/utils.cpp src/PhasedArrayGroundStation.cpp)
add_executable(utils_tests tests/utils_tests.cpp src/Satellite.cpp src/utils.cpp src/PhasedArrayGroundStation.cpp)
add_executable(gs_tests tests/PhasedArrayGroundStation_tests.cpp src/Satellite.cpp src/utils.cpp src/PhasedArrayGroundStation.cpp)

target_link_libraries(run PRIVATE nlohmann_json::nlohmann_json Eigen3::Eigen)
target_link_libraries(Satellite_tests PRIVATE nlohmann_json::nlohmann_json gtest_main Eigen3::Eigen)
target_link_libraries(utils_tests PRIVATE nlohmann_json::nlohmann_json gtest_main Eigen3::Eigen)
target_link_libraries(gs_tests PRIVATE nlohmann_json::nlohmann_json gtest_main Eigen3::Eigen)
