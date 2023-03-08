include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

macro(build_tests testname)
    add_executable(${testname} ${ARGN})
    link_library(${testname} fftpack ${PROJECT_INCLUDE_DIR})
    link_library(${testname} test-drive ${test-drive_INCLUDE_DIR})
    add_test(
        NAME ${testname}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:${testname}>
    )
endmacro()

set(FFTPACK_TEST_SOURCES
    test_fftpack_dct.f90
    test_fftpack_fft.f90
    test_fftpack_rfft.f90
    test_fftpack_utils.f90
    test_fftpack.f90
)

# Run the original FFTPACK test
build_tests(fftpack_original_tests tstfft.f)

# Run the additional FFTPACK tests
build_tests(fftpack_tests ${FFTPACK_TEST_SOURCES})
