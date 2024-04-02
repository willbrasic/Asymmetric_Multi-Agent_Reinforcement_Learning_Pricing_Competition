cmake_minimum_required(VERSION 3.10)

project(SARSA_Qlearning_Heat_Map)  # Name of your project

# Add the executable target
add_executable(SARSA_Qlearning_Heat_Map main.cpp Functions.cpp Global.cpp)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Get compiler to optimize the most
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")