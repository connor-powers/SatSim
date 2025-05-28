cmake_minimum_required(VERSION 3.10)
project(example_satsim_project)
include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4
)
FetchContent_MakeAvailable(json Eigen)

find_library(SATSIM NAMES SatSim)

message(STATUS ${CMAKE_INSTALL_PREFIX}/lib)

add_executable(run simulation_setup.cpp)
add_executable(sim_profiling_executable sim_profiling.cpp)
target_include_directories(run PUBLIC /usr/local/include/SatSim)
target_include_directories(sim_profiling_executable PUBLIC /usr/local/include/SatSim)
target_link_libraries(run PRIVATE nlohmann_json::nlohmann_json Eigen3::Eigen ${SATSIM})
target_link_libraries(sim_profiling_executable PRIVATE nlohmann_json::nlohmann_json Eigen3::Eigen ${SATSIM})
set_target_properties(run PROPERTIES
    BUILD_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"    
    )
set_target_properties(sim_profiling_executable PROPERTIES
    BUILD_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"    
    )