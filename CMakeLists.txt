cmake_minimum_required(VERSION 3.10.0)
project(mpc_solver VERSION 0.1.0 LANGUAGES C CXX)

include_directories(include)

set(source_files
    src/main.cpp
    src/dynamic_game_planner.cpp
    src/utils.cpp
    src/vehicle_state.cpp
)

add_executable(mpc_solver ${source_files})

