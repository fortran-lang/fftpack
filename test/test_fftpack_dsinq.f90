program tester

    call test_fftpack_dsinq_real
    print *, "All tests in `test_fftpack_dsinq` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dsinq_real
        use fftpack, only: dsinqi, dsinqf, dsinqb
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: w(3*4 + 15)
        real(kind=dp) :: x(4) = [1, 2, 3, 4]
        real(kind=dp) :: eps = 1.0e-10_dp

        call dsinqi(4, w)
        call dsinqf(4, x, w)
        call check(sum(abs(x - [13.137071184544091_dp, -1.6199144044217753_dp, &
                                0.72323134608584416_dp, -0.51978306494828974_dp])) < eps, msg="`dsinqf` failed.")

        call dsinqb(4, x, w)
        call check(sum(abs(x/(4.0_dp*4.0_dp) - [1.0000000000000002_dp, 2.0_dp, 3.0000000000000009_dp, 4.0_dp])) < eps, &
                   msg="`dsinqb` failed.")

    end subroutine test_fftpack_dsinq_real

end program tester
