add_executable(bench1 bench01_zfft.f90)
target_link_libraries(bench1 fftpack)

add_executable(bench2 bench02_zfft.f90)
target_link_libraries(bench2 fftpack)

add_executable(bench3 bench03_dfft.f90)
target_link_libraries(bench3 fftpack)

add_executable(rfft_example)
target_link_libraries(rfft_example fftpack)

add_executable(complex_transforms)
target_link_libraries(complex_transforms fftpack)

add_executable(real_transforms)
target_link_libraries(real_transforms fftpack)
