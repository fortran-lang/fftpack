if (NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

if(WIN32)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    add_link_options(/Qoption,link,/STACK:8388608)
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    add_link_options(-Wl,--stack,8388608)
  endif()
endif()

macro(ADDTEST name)
    add_executable(test_${name} test_${name}.f90)
    target_link_libraries(test_${name} "${PROJECT_NAME}" "test-drive::test-drive")
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:test_${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDTEST)

macro(ADDTESTPP name)
    add_executable(test_${name} test_${name}.F90)
    target_link_libraries(test_${name} "${PROJECT_NAME}" "test-drive::test-drive")
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:test_${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDTESTPP)

ADDTEST(fftpack_fft)
ADDTEST(fftpack_rfft)
ADDTEST(fftpack_dct)
ADDTEST(fftpack_utils)
ADDTEST(fftpack_original)
