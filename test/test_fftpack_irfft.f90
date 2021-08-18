program tester

    call test_fftpack_irfft()
    print *, "All tests in `test_fftpack_irfft` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_irfft
        use fftpack, only: rfft, irfft
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10

        real(kind=dp) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(irfft(rfft(x))/4.0 - [real(kind=dp) :: 1, 2, 3, 4])) < eps, &
                   msg="sum(abs(irfft(rfft(x))/4.0 - [real(kind=dp) :: 1, 2, 3, 4])) < eps failed.")
        call check(sum(abs(irfft(rfft(x), 2) - [real(kind=dp) :: 8, 12])) < eps, &
                   msg="sum(abs(irfft(rfft(x), 2) - [real(kind=dp) :: 8, 12])) < eps failed.")
        call check(sum(abs(irfft(rfft(x, 2), 4) - [real(kind=dp) :: 1, 3, 5, 3])) < eps, &
                   msg="sum(abs(irfft(rfft(x, 2), 4) - [real(kind=dp) :: 1, 3, 5, 3])) failed.")

    end subroutine test_fftpack_irfft

end program tester
