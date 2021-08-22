program tester

    call test_fftpack_iqct()
    print *, "All tests in `test_fftpack_iqct` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_iqct
        use fftpack, only: qct, iqct
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10

        real(kind=dp) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(iqct(qct(x))/(4.0_dp*4.0_dp) - [real(kind=dp) :: 1, 2, 3, 4])) < eps, &
                   msg="`iqct(qct(x)/(4.0_dp*4.0_dp)` failed.")
        call check(sum(abs(iqct(qct(x), 2)/(4.0_dp*2.0_dp) - [1.4483415291679655_dp, 7.4608849947753271_dp])) < eps, &
                   msg="`iqct(qct(x), 2)/(4.0_dp*2.0_dp)` failed.")
        call check(sum(abs(iqct(qct(x, 2), 4)/(4.0_dp*4.0_dp) - [0.5_dp, 0.70932417358418376_dp, 1.0_dp, &
                                                                 0.78858050747473762_dp])) < eps, &
                   msg="`iqct(qct(x, 2), 4)/(4.0_dp*4.0_dp)` failed.")

    end subroutine test_fftpack_iqct

end program tester
