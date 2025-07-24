macro(ADD_EXAMPLE name)
    add_executable(example_${name} example_${name}.f90)
    target_link_libraries(example_${name} "${PROJECT_NAME}")
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:example_${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADD_EXAMPLE)

ADD_EXAMPLE(bench01_zfft)
ADD_EXAMPLE(bench02_zfft)
ADD_EXAMPLE(bench03_dfft)
ADD_EXAMPLE(complex_transforms)
ADD_EXAMPLE(real_transforms)
