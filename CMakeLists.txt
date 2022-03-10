cmake_minimum_required(VERSION 3.20)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

project(Polynomial)

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories( $(Polynomial) )

add_executable(Polynomial/polynomial.h Polynomial/polynomial.cpp)

install(DIRECTORY Polynomial/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Polynomial)
