# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# By default, static library
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Export all symbols on Windows when building libraries
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Utilize the GNU installation structure
include(GNUInstallDirs)

# Locate the local include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include)
set(PROJECT_INCLUDE_DIR ${PROJECT_INCLUDE_DIR} PARENT_SCOPE)