include(FetchContent)

# If found, use the installed version; else, import the library
if (${test-drive_FOUND})
    # Inform the user of what's going on
    message(STATUS "TEST-DRIVE (${test-drive_VERSION}) library found.")

    # Get the mod file location
    get_target_property(test-drive_INCLUDE_DIR test-drive INTERFACE_INCLUDE_DIRECTORIES)
else()
    # Inform the user of what's going on
    message(STATUS "TEST-DRIVE library not found.  Downloading appropriate repository.")
    
    # Fetch the proper content
    FetchContent_Declare(
        test-drive
        GIT_TAG "origin/main"
        GIT_REPOSITORY "https://github.com/fortran-lang/test-drive"
        OVERRIDE_FIND_PACKAGE
    )

    FetchContent_MakeAvailable(test-drive)
    set(test-drive_INCLUDE_DIR ${test-drive_BINARY_DIR}/include)
endif()

# Make a parent-scope variable for the library
set(test-drive_LIBRARY test-drive)
set(test-drive_LIBRARY ${test-drive_LIBRARY} PARENT_SCOPE)

# Make a parent-scope variable locating the include directory for test-drive
set(test-drive_INCLUDE_DIR ${test-drive_INCLUDE_DIR} PARENT_SCOPE)
