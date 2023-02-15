# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(CMakePackageConfigHelpers)

# Install the library and necessary include files
install_library(
    ${PROJECT_NAME}
    ${CMAKE_INSTALL_LIBDIR}
    ${CMAKE_INSTALL_BINDIR}
    ${PROJECT_INCLUDE_DIR}
    ${CMAKE_INSTALL_PREFIX}
)

# Define the version file
write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(
    EXPORT ${PROJECT_NAME}Targets
    FILE ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Targets.cmake
)

# Define the configuration file
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    COPYONLY
)

install(
    EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
install(
    FILES
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/template.pc
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc
  @ONLY
)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)