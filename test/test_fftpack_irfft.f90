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
        use fftpack, only: rfft, irfft, dp
        real(kind=dp) :: eps = 1.0e-10_dp

        real(kind=dp) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(irfft(rfft(x))/4.0_dp - [real(kind=dp) :: 1, 2, 3, 4])) < eps, &
                   msg="`irfft(rfft(x))/4.0_dp` failed.")
        call check(sum(abs(irfft(rfft(x), 2) - [real(kind=dp) :: 8, 12])) < eps, &
                   msg="`irfft(rfft(x), 2)` failed.")
        call check(sum(abs(irfft(rfft(x, 2), 4) - [real(kind=dp) :: 1, 3, 5, 3])) < eps, &
                   msg="`irfft(rfft(x, 2), 4)` failed.")

    end subroutine test_fftpack_irfft

end program tester
