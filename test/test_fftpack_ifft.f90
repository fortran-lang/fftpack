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
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10_dp

        complex(kind=dp) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(ifft(fft(x))/4.0_dp - [complex(kind=dp) :: 1, 2, 3, 4])) < eps, &
                   msg="`ifft(fft(x))/4.0_dp` failed.")
        call check(sum(abs(ifft(fft(x), 2) - [complex(kind=dp) ::(8, 2), (12, -2)])) < eps, &
                   msg="`ifft(fft(x), 2)` failed.")
        call check(sum(abs(ifft(fft(x, 2), 4) - [complex(kind=dp) ::(2, 0), (3, -1), (4, 0), (3, 1)])) < eps, &
                   msg="`ifft(fft(x, 2), 4)` failed.")

    end subroutine test_fftpack_ifft

end program tester
