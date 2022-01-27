program tester

    call test_fftpack_dcost_real
    print *, "All tests in `test_fftpack_dcost` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dcost_real
        use fftpack, only: dcosti, dcost
        use fftpack_kind
        real(kind=rk) :: w(3*4 + 15)
        real(kind=rk) :: x(4) = [1, 2, 3, 4]
        real(kind=rk) :: eps = 1.0e-10_rk

        call dcosti(4, w)
        call dcost(4, x, w)
        call check(all(x == [real(kind=rk) :: 15, -4, 0, -1.0000000000000009_rk]), msg="`dcosti` failed.")

        call dcost(4, x, w)
        call check(all(x/(2.0_rk*(4.0_rk - 1.0_rk)) == [real(kind=rk) :: 1, 2, 3, 4]), msg="`dcost` failed.")

    end subroutine test_fftpack_dcost_real

end program tester
