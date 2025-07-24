macro(ADD_EXAMPLE name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} "${PROJECT_NAME}")
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADD_EXAMPLE)

ADD_EXAMPLE(bench01_zfft)
ADD_EXAMPLE(bench02_zfft)
ADD_EXAMPLE(bench03_dfft)
ADD_EXAMPLE(complex_transforms)
ADD_EXAMPLE(real_transforms)
