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
      use fftpack_kind

    use fftpack_kind
        use fftpack, only: zffti, zfftf, zfftb
    use fftpack_kind

        complex(kind=rk) :: x(4) = [1, 2, 3, 4]
        real(kind=rk) :: w(31)

        call zffti(4, w)
        call zfftf(4, x, w)
        call check(all(x == [complex(kind=rk) ::(10, 0), (-2, 2), (-2, 0), (-2, -2)]), &
                   msg="`zfftf` failed.")

        call zfftb(4, x, w)
        call check(all(x/4.0_rk == [complex(kind=rk) ::(1, 0), (2, 0), (3, 0), (4, 0)]), &
                   msg="`zfftb` failed.")

    end subroutine test_fftpack_zfft

end program tester
