program tester

    call test_fftpack_fft()
    print *, "All tests in `test_fftpack_fft` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_fft
        use fftpack, only: fft
    use fftpack_kind
        real(kind=rk) :: eps = 1.0e-10_rk

        complex(kind=rk) :: x(3) = [1.0_rk, 2.0_rk, 3.0_rk]

        call check(sum(abs(fft(x, 2) - [(3.0_rk, 0.0_rk), (-1.0_rk, 0.0_rk)])) < eps, &
                   msg="`fft(x, 2)` failed.")
        call check(sum(abs(fft(x, 3) - fft(x))) < eps, &
                   msg="`fft(x, 3)` failed.")
        call check(sum(abs(fft(x, 4) - [(6.0_rk, 0.0_rk), (-2.0_rk, -2.0_rk), (2.0_rk, 0.0_rk), (-2.0_rk, 2.0_rk)])) < eps, &
                   msg="`fft(x, 4)` failed.")

    end subroutine test_fftpack_fft

end program tester
