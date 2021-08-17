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
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10

        complex(kind=dp) :: x(3)

        x = [real(kind=dp) :: 1.0, 2.0, 3.0]

        call check(sum(abs(fft(x, 2) - [complex(kind=dp) ::(3.0, 0.0), (-1.0, 0.0)])) < eps, &
                   msg="sum(abs(fft(x,2) - [complex(kind=dp) ::(3.0, 0.0), (-1.0, 0.0)])) < eps failed.")
        call check(sum(abs(fft(x, 3) - fft(x))) < eps, &
                   msg="sum(abs(fft(x,3) - fft(3))) < eps failed.")
        call check(sum(abs(fft(x, 4) - [complex(kind=dp) ::(6.0, 0.0), (-2.0, -2.0), (2.0, 0.0), (-2.0, 2.0)])) < eps, &
                   msg="sum(abs(fft(x,4) - [complex(kind=dp) ::(6.0, 0.0), (-2.0, -2.0), (2.0,0.0), (-2.0,2.0)])) < eps failed.")

    end subroutine test_fftpack_fft

end program tester
