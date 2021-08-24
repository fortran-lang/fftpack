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
        use iso_fortran_env, only: dp => real64

        real(kind=dp) :: x(3) = [9, -9, 3]

        call check(all(dct(x, 2) == [real(kind=dp) :: 0, 18]), msg="`dct(x, 2)` failed.")
        call check(all(dct(x, 3) == dct(x)), msg="`dct(x, 3)` failed.")
        call check(all(dct(x, 4) == [real(kind=dp) :: -3, -3.0000000000000036_dp, 15, 33]), msg="`dct(x, 4)` failed.")

    end subroutine test_fftpack_dct

    subroutine test_fftpack_idct
        use fftpack, only: dct, idct
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10_dp
        real(kind=dp) :: x(4) = [1, 2, 3, 4]

        call check(all(idct(dct(x))/(2.0_dp*(4.0_dp - 1.0_dp)) == [real(kind=dp) :: 1, 2, 3, 4]), &
                   msg="`idct(dct(x))/(2.0_dp*(4.0_dp-1.0_dp))` failed.")
        call check(all(idct(dct(x), 2)/(2.0_dp*(2.0_dp - 1.0_dp)) == [real(kind=dp) :: 5.5, 9.5]), &
                   msg="`idct(dct(x), 2)/(2.0_dp*(2.0_dp-1.0_dp))` failed.")
        call check(all(idct(dct(x, 2), 4)/(2.0_dp*(4.0_dp - 1.0_dp)) == &
                   [0.16666666666666666_dp, 0.33333333333333331_dp, 0.66666666666666663_dp, 0.83333333333333315_dp]), &
                   msg="`idct(dct(x, 2), 4)/(2.0_dp*(4.0_dp-1.0_dp))` failed.")

    end subroutine test_fftpack_idct

end program tester
