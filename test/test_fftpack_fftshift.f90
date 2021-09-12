program tester

    call test_fftpack_fftshift_complex
    call test_fftpack_fftshift_real
    print *, "All tests in `test_fftpack_fftshift` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_fftshift_complex
        use fftpack, only: fftshift
    use fftpack_kind

        complex(kind=rk) :: xeven(4) = [1, 2, 3, 4]
        complex(kind=rk) :: xodd(5) = [1, 2, 3, 4, 5]

        call check(all(fftshift(xeven) == [complex(kind=rk) :: 3, 4, 1, 2]), &
                   msg="all(fftshift(xeven) == [complex(kind=rk) :: 3, 4, 1, 2]) failed.")
        call check(all(fftshift(xodd) == [complex(kind=rk) :: 4, 5, 1, 2, 3]), &
                   msg="all(fftshift(xodd) == [complex(kind=rk) :: 4, 5, 1, 2, 3]) failed.")

    end subroutine test_fftpack_fftshift_complex

    subroutine test_fftpack_fftshift_real
        use fftpack, only: fftshift
    use fftpack_kind

        real(kind=rk) :: xeven(4) = [1, 2, 3, 4]
        real(kind=rk) :: xodd(5) = [1, 2, 3, 4, 5]

        call check(all(fftshift(xeven) == [real(kind=rk) :: 3, 4, 1, 2]), &
                   msg="all(fftshift(xeven) == [real(kind=rk) :: 3, 4, 1, 2]) failed.")
        call check(all(fftshift(xodd) == [real(kind=rk) :: 4, 5, 1, 2, 3]), &
                   msg="all(fftshift(xodd) == [real(kind=rk) :: 4, 5, 1, 2, 3]) failed.")

    end subroutine test_fftpack_fftshift_real

end program tester
