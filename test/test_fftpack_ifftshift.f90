program tester

    call test_fftpack_ifftshift_complex
    call test_fftpack_ifftshift_real
    print *, "All tests in `test_fftpack_ifftshift` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_ifftshift_complex
        use fftpack, only: ifftshift
    use fftpack_kind
        integer :: i

        complex(kind=rk) :: xeven(4) = [3, 4, 1, 2]
        complex(kind=rk) :: xodd(5) = [4, 5, 1, 2, 3]

        call check(all(ifftshift(xeven) == [complex(kind=rk) ::(i, i=1, 4)]), &
                   msg="all(ifftshift(xeven) == [complex(kind=rk) ::(i, i=1, 4)]) failed.")
        call check(all(ifftshift(xodd) == [complex(kind=rk) ::(i, i=1, 5)]), &
                   msg="all(ifftshift(xodd) == [complex(kind=rk) ::(i, i=1, 5)]) failed.")

    end subroutine test_fftpack_ifftshift_complex

    subroutine test_fftpack_ifftshift_real
        use fftpack, only: ifftshift
    use fftpack_kind
        integer :: i

        real(kind=rk) :: xeven(4) = [3, 4, 1, 2]
        real(kind=rk) :: xodd(5) = [4, 5, 1, 2, 3]

        call check(all(ifftshift(xeven) == [real(kind=rk) ::(i, i=1, 4)]), &
                   msg="all(ifftshift(xeven) == [real(kind=rk) ::(i, i=1, 4)]) failed.")
        call check(all(ifftshift(xodd) == [real(kind=rk) ::(i, i=1, 5)]), &
                   msg="all(ifftshift(xodd) == [real(kind=rk) ::(i, i=1, 5)]) failed.")

    end subroutine test_fftpack_ifftshift_real

end program tester
