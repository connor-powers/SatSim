cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
include_directories(include)
include_directories(/usr/local/eigen-3.4.0) #change this based on where your Eigen installation is

project(orbital_sim)
include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

add_executable(run src/orbits_core.cpp src/utils.cpp src/Satellite.cpp)
target_link_libraries(run PRIVATE nlohmann_json::nlohmann_json)