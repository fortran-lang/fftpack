program tester

    call test_fftpack_qst()
    print *, "All tests in `test_fftpack_qst` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_qst
        use fftpack, only: qst
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10_dp

        real(kind=dp) :: x(3) = [9, -9, 3]

        call check(sum(abs(qst(x, 2) - [3.7279220613578570_dp, 21.727922061357859_dp])) < eps, &
                   msg="`qst(x, 2)` failed.")
        call check(sum(abs(qst(x, 3) - qst(x))) < eps, &
                   msg="`qst(x,3)` failed.")
        call check(sum(abs(qst(x, 4) - [-0.29634308371852036_dp, 1.6058089296547653_dp, &
                                        27.061653052370481_dp, 25.159501038997192_dp])) < eps, &
                   msg="`qst(x, 4)` failed.")

    end subroutine test_fftpack_qst

end program tester
