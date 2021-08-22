program tester

    call test_fftpack_dcosq_real
    print *, "All tests in `test_fftpack_dcosq` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dcosq_real
        use fftpack, only: dcosqi, dcosqf, dcosqb
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: w(3*4 + 15)
        real(kind=dp) :: x(4) = [1, 2, 3, 4]
        real(kind=dp) :: eps = 1.0e-10_dp

        call dcosqi(4, w)
        call dcosqf(4, x, w)
        call check(sum(abs(x - [11.999626276085150_dp, -9.1029432177492193_dp, &
                                2.6176618435106480_dp, -1.5143449018465791_dp])) < eps, msg="`dcosqf` failed.")
        call dcosqb(4, x, w)    !!
        call check(sum(abs(x/(4.0_dp*4.0_dp) - [real(kind=dp) :: 1, 2, 3, 4])) < eps, msg="`dcosqb` failed.")

    end subroutine test_fftpack_dcosq_real

end program tester
