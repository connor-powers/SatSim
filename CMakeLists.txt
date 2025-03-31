cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
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



add_executable(run src/simulation_setup.cpp src/utils.cpp src/Satellite.cpp)
add_executable(circular_orbit_tests tests/circular_orbit_tests.cpp src/Satellite.cpp src/utils.cpp)
add_executable(elliptical_orbit_tests tests/elliptical_orbit_tests.cpp src/Satellite.cpp src/utils.cpp)
add_executable(attitude_tests tests/attitude_tests.cpp src/Satellite.cpp src/utils.cpp)

target_link_libraries(run PRIVATE nlohmann_json::nlohmann_json Eigen3::Eigen)
target_link_libraries(circular_orbit_tests PRIVATE nlohmann_json::nlohmann_json gtest_main Eigen3::Eigen)
target_link_libraries(elliptical_orbit_tests PRIVATE nlohmann_json::nlohmann_json gtest_main Eigen3::Eigen)
target_link_libraries(attitude_tests PRIVATE nlohmann_json::nlohmann_json gtest_main Eigen3::Eigen)
