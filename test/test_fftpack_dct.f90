program tester

    call test_fftpack_dct()
    call test_fftpack_idct()
    print *, "All tests in `test_fftpack_dct` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dct
        use fftpack, only: dct
        use fftpack_kind

        real(kind=rk) :: x(3) = [9, -9, 3]

        call check(all(dct(x, 2) == [real(kind=rk) :: 0, 18]), msg="`dct(x, 2)` failed.")
        call check(all(dct(x, 3) == dct(x)), msg="`dct(x, 3)` failed.")
        call check(all(dct(x, 4) == [real(kind=rk) :: -3, -3.0000000000000036_rk, 15, 33]), msg="`dct(x, 4)` failed.")

    end subroutine test_fftpack_dct

    subroutine test_fftpack_idct
        use fftpack, only: dct, idct
        use iso_fortran_env, only: rk => real64
        real(kind=rk) :: eps = 1.0e-10_rk
        real(kind=rk) :: x(4) = [1, 2, 3, 4]

        call check(all(idct(dct(x))/(2.0_rk*(4.0_rk - 1.0_rk)) == [real(kind=rk) :: 1, 2, 3, 4]), &
                   msg="`idct(dct(x))/(2.0_rk*(4.0_rk-1.0_rk))` failed.")
        call check(all(idct(dct(x), 2)/(2.0_rk*(2.0_rk - 1.0_rk)) == [real(kind=rk) :: 5.5, 9.5]), &
                   msg="`idct(dct(x), 2)/(2.0_rk*(2.0_rk-1.0_rk))` failed.")
        call check(all(idct(dct(x, 2), 4)/(2.0_rk*(4.0_rk - 1.0_rk)) == &
                   [0.16666666666666666_rk, 0.33333333333333331_rk, 0.66666666666666663_rk, 0.83333333333333315_rk]), &
                   msg="`idct(dct(x, 2), 4)/(2.0_rk*(4.0_rk-1.0_rk))` failed.")

    end subroutine test_fftpack_idct

end program tester
