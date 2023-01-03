# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
  fftpack
  LANGUAGES Fortran
  VERSION 1.0.0
)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Confgiure everything
add_subdirectory(configure)

# Source
add_subdirectory(src)
add_fortran_library(
  ${PROJECT_NAME}
  ${PROJECT_INCLUDE_DIR}
  ${CMAKE_INSTALL_INCLUDEDIR}
  ${PROJECT_VERSION}
  ${PROJECT_VERSION_MAJOR}
  ${FFTPACK_SOURCES}
)

# Installation
add_subdirectory(install)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build FFTPACK tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
  enable_testing()
  add_subdirectory(dependencies)
  add_subdirectory(test)
endif()

# Examples
option(BUILD_FFTPACK_EXAMPLES "Build FFTPACK examples")
message(STATUS "Build FFTPACK examples: ${BUILD_FFTPACK_EXAMPLES}")
if (BUILD_FFTPACK_EXAMPLES)
  add_subdirectory(example)
endif()