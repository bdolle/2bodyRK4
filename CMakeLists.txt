cmake_minimum_required(VERSION 3.9)
project(twobodytest)

set(CMAKE_CXX_STANDARD 11)

add_executable(twobodytest main.cpp Coords.cpp Coords.h Body.cpp Body.h SolarSystem.cpp SolarSystem.h)