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
        real(kind=dp) :: eps = 1.0e-10

        complex(kind=dp) :: x(4)

        x = [real(kind=dp) :: 1.0, 2.0, 3.0, 4.0]
        
        call check(sum(abs(ifft(fft(x))/4.0 - [complex(kind=dp) ::(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)])) < eps, &
        msg="abs(sum(ifft(fft(x))/4.0 - [complex(kind=dp) ::(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)])) < eps failed.")
        call check(sum(abs(ifft(fft(x),2) - [complex(kind=dp) ::(8.0, 2.0), (12.0, -2.0)])) < eps, &
                   msg="abs(sum(ifft(fft(x),2) - [complex(kind=dp) ::(8.0, 2.0), (12.0, -2.0)])) < eps failed.")
        call check(sum(abs(ifft(fft(x,2),4) - [complex(kind=dp) ::(2.0, 0.0), (3.0, -1.0), (4.0,0.0), (3.0,1.0)])) < eps, &
        msg="abs(sum(ifft(fft(x,2),4) - [complex(kind=dp) ::(2.0, 0.0), (3.0, -1.0), (4.0,0.0), (3.0,1.0)])) < eps failed.")

    end subroutine test_fftpack_ifft

end program tester
