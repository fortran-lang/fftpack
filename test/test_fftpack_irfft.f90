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
    use fftpack_kind
        real(kind=rk) :: eps = 1.0e-10_rk

        real(kind=rk) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(irfft(rfft(x))/4.0_rk - [real(kind=rk) :: 1, 2, 3, 4])) < eps, &
                   msg="`irfft(rfft(x))/4.0_rk` failed.")
        call check(sum(abs(irfft(rfft(x), 2) - [real(kind=rk) :: 8, 12])) < eps, &
                   msg="`irfft(rfft(x), 2)` failed.")
        call check(sum(abs(irfft(rfft(x, 2), 4) - [real(kind=rk) :: 1, 3, 5, 3])) < eps, &
                   msg="`irfft(rfft(x, 2), 4)` failed.")

    end subroutine test_fftpack_irfft

end program tester
