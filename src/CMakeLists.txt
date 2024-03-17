# Locate the source directories
set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(subdir "${dir}/fftpack")

# The source files
set(FFTPACK_SOURCES
    ${subdir}/cfftb1.f90
    ${subdir}/cfftf1.f90
    ${subdir}/cffti1.f90
    ${subdir}/cosqb1.f90
    ${subdir}/cosqf1.f90
    ${subdir}/dcosqb.f90
    ${subdir}/dcosqf.f90
    ${subdir}/dcosqi.f90
    ${subdir}/dcost.f90
    ${subdir}/dcosti.f90
    ${subdir}/dfftb.f90
    ${subdir}/dfftf.f90
    ${subdir}/dffti.f90
    ${subdir}/dsinqb.f90
    ${subdir}/dsinqf.f90
    ${subdir}/dsinqi.f90
    ${subdir}/dsint.f90
    ${subdir}/dsinti.f90
    ${subdir}/dzfftb.f90
    ${subdir}/dzfftf.f90
    ${subdir}/dzffti.f90
    ${subdir}/ezfft1.f90
    ${dir}/fftpack.f90
    ${subdir}/fftpack_dct.f90
    ${subdir}/fftpack_fft.f90
    ${subdir}/fftpack_fftshift.f90
    ${subdir}/fftpack_ifft.f90
    ${subdir}/fftpack_ifftshift.f90
    ${subdir}/fftpack_irfft.f90
    ${subdir}/fftpack_rfft.f90
    ${subdir}/fftpack_utils.f90
    ${subdir}/passb.f90
    ${subdir}/passb2.f90
    ${subdir}/passb3.f90
    ${subdir}/passb4.f90
    ${subdir}/passb5.f90
    ${subdir}/passf.f90
    ${subdir}/passf2.f90
    ${subdir}/passf3.f90
    ${subdir}/passf4.f90
    ${subdir}/passf5.f90
    ${subdir}/radb2.f90
    ${subdir}/radb3.f90
    ${subdir}/radb4.f90
    ${subdir}/radb5.f90
    ${subdir}/radbg.f90
    ${subdir}/radf2.f90
    ${subdir}/radf3.f90
    ${subdir}/radf4.f90
    ${subdir}/radf5.f90
    ${subdir}/radfg.f90
    ${subdir}/rfftb1.f90
    ${subdir}/rfftf1.f90
    ${subdir}/rffti1.f90
    ${subdir}/rk.f90
    ${subdir}/sint1.f90
    ${subdir}/zfftb.f90
    ${subdir}/zfftf.f90
    ${subdir}/zffti.f90
)
set(FFTPACK_SOURCES ${FFTPACK_SOURCES} PARENT_SCOPE)
