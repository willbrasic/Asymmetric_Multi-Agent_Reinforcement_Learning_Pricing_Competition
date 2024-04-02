cmake_minimum_required(VERSION 3.0)  # Minimum version of CMake required

project(SARSA_Qlearning)  # Name of your project

# Add the executable target
add_executable(SARSA_Qlearning main.cpp Functions.cpp Global.cpp Agent.cpp SARSA_Agent.cpp Qlearning_Agent.cpp)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Get compiler to optimize the most
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
