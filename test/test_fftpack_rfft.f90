program tester

    call test_fftpack_rfft()
    print *, "All tests in `test_fftpack_rfft` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_rfft
        use fftpack, only: rfft
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10

        real(kind=dp) :: x(3) = [9, -9, 3]

        call check(sum(abs(rfft(x, 2) - [real(kind=dp) :: 0, 18])) < eps, &
                   msg="sum(abs(rfft(x,2) - [real(kind=dp) :: 0, 18])) < eps failed.")
        call check(sum(abs(rfft(x, 3) - rfft(x))) < eps, &
                   msg="sum(abs(rfft(x,3) - rfft(3))) < eps failed.")
        call check(sum(abs(rfft(x, 4) - [real(kind=dp) :: 3, 6, 9, 21])) < eps, &
                   msg="sum(abs(rfft(x,4) - [real(kind=dp) ::3, 6, 9, 21])) < eps failed.")

    end subroutine test_fftpack_rfft

end program tester
