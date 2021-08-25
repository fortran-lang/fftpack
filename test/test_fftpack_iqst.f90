program tester

    call test_fftpack_iqst()
    print *, "All tests in `test_fftpack_iqst` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_iqst
        use fftpack, only: qst, iqst
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10_dp

        real(kind=dp) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(iqst(qst(x))/(4.0_dp*4.0_dp) - [1.0000000000000002_dp, 2.0_dp, &
                                                           3.0000000000000009_dp, 4.0_dp])) < eps, &
                   msg="`iqst(qst(x)/(4.0_dp*4.0_dp)` failed.")
        call check(sum(abs(iqst(qst(x), 2)/(4.0_dp*2.0_dp) - [4.0719298296065567_dp, 7.3784927944829333_dp])) < eps, &
                   msg="`iqst(qst(x), 2)/(4.0_dp*2.0_dp)` failed.")
        call check(sum(abs(iqst(qst(x, 2), 4)/(4.0_dp*4.0_dp) - [0.19134171618254481_dp, 0.5_dp, &
                                                                 0.84462319862073310_dp, 1.0_dp])) < eps, &
                   msg="`iqst(qst(x, 2), 4)/(4.0_dp*4.0_dp)` failed.")

    end subroutine test_fftpack_iqst

end program tester
