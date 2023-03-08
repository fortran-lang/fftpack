program test_fftpack
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type
    use test_fftpack_fft, only: collect_fft
    use test_fftpack_rfft, only: collect_rfft
    use test_fftpack_dct, only: collect_dct
    use test_fftpack_utils, only: collect_utils
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
                 new_testsuite("fft", collect_fft), &
                 new_testsuite("rfft", collect_rfft), &
                 new_testsuite("dct", collect_dct), &
                 new_testsuite("utils", collect_utils) &
                 ]

    do is = 1, size(testsuites)
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

end program test_fftpack
