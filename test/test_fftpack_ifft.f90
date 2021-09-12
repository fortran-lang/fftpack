program tester

    call test_fftpack_ifft()
    print *, "All tests in `test_fftpack_ifft` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_ifft
        use fftpack, only: fft, ifft
    use fftpack_kind
        real(kind=rk) :: eps = 1.0e-10_rk

        complex(kind=rk) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(ifft(fft(x))/4.0_rk - [complex(kind=rk) :: 1, 2, 3, 4])) < eps, &
                   msg="`ifft(fft(x))/4.0_rk` failed.")
        call check(sum(abs(ifft(fft(x), 2) - [complex(kind=rk) ::(8, 2), (12, -2)])) < eps, &
                   msg="`ifft(fft(x), 2)` failed.")
        call check(sum(abs(ifft(fft(x, 2), 4) - [complex(kind=rk) ::(2, 0), (3, -1), (4, 0), (3, 1)])) < eps, &
                   msg="`ifft(fft(x, 2), 4)` failed.")

    end subroutine test_fftpack_ifft

end program tester
