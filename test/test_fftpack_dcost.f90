program tester

    call test_fftpack_dcost_real
    print *, "All tests in `test_fftpack_dcost` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dcost_real
        use fftpack, only: dcosti, dcost
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: w(3*4 + 15)
        real(kind=dp) :: x(4) = [1, 2, 3, 4]
        real(kind=dp) :: eps = 1.0e-10_dp

        call dcosti(4, w)
        call dcost(4, x, w)
        call check(all(x == [real(kind=dp) :: 15, -4, 0, -1.0000000000000009_dp]), msg="`dcosti` failed.")

        call dcost(4, x, w)
        call check(all(x/(2.0_dp*(4.0_dp - 1.0_dp)) == [real(kind=dp) :: 1, 2, 3, 4]), msg="`dcost` failed.")

    end subroutine test_fftpack_dcost_real

end program tester
