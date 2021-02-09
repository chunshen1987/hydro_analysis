cmake_minimum_required (VERSION 3.0)
project (hydro_analysis)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")
set (CMAKE_CXX_COMPILER "h5c++")
set (CMAKE_CXX_FLAGS "-g -O3 -Wall -std=c++11")

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

add_subdirectory (src)
