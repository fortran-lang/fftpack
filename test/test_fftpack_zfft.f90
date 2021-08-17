program tester

    call test_fftpack_zfft()
    print *, "All tests in `test_fftpack_zfft` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_zfft()
        use fftpack, only: zffti, zfftf, zfftb
        use iso_fortran_env, only: dp => real64

        complex(kind=dp) :: x(4)
        real(kind=dp) :: w(31)

        x = [real(kind=dp) :: 1.0, 2.0, 3.0, 4.0]
        call zffti(4, w)
        call zfftf(4, x, w)
        call check(all(x == [complex(kind=dp) ::(10.0, 0.0), (-2.0, 2.0), (-2.0, 0.0), (-2.0, -2.0)]), &
                   msg="all(x==[complex(kind=dp) :: (10.0,0.0), (-2.0,2.0), (-2.0,0.0), (-2.0,-2.0)]) failed.")

        call zfftb(4, x, w)
        call check(all(x/4.0_dp == [complex(kind=dp) ::(1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)]), &
                   msg="all(x/4==[complex(kind=dp) :: (1.0,0.0), (2.0,0.0), (3.0,0.0), (4.0,0.0)]) failed.")

    end subroutine test_fftpack_zfft

end program tester
