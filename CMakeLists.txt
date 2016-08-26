cmake_minimum_required (VERSION 2.6)
project (hydro_analysis)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")
set (CMAKE_CXX_COMPILER "h5c++")
set (CMAKE_CXX_FLAGS "-O3 -Wall")

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

add_subdirectory (src)
