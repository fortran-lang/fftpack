program tester

    call test_fftpack_dst()
    call test_fftpack_idst()
    print *, "All tests in `test_fftpack_dst` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dst
        use fftpack, only: dst
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10_dp
        real(kind=dp) :: x(3) = [1, 2, 3]

        call check(sum(abs(dst(x, 2) - [5.1961524227066320_dp, -1.7320508075688772_dp])) < eps, &
                   msg="`dst(x, 2)` failed.")
        call check(all(dst(x, 3) == dst(x)), msg="`dst(x, 3)` failed.")
        call check(sum(abs(dst(x, 4) - [10.686135667536481_dp, 0.72654252800536079_dp, &
                                        -3.9757394903344245_dp, 3.0776835371752531_dp])) < eps, &
                   msg="`dst(x, 4)` failed.")

    end subroutine test_fftpack_dst

    subroutine test_fftpack_idst
        use fftpack, only: dst, idst
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: eps = 1.0e-10_dp
        real(kind=dp) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(idst(dst(x))/(2.0_dp*(4.0_dp + 1.0_dp)) - [real(kind=dp) :: 1, 2, 3, 4])) < eps, &
                   msg="`idst(dst(x))/(2.0_dp*(4.0_dp+1.0_dp))` failed.")
        call check(sum(abs(idst(dst(x), 2)/(2.0_dp*(2.0_dp + 1.0_dp)) - [2.4556173659421150_dp, 6.4288897274009438_dp])) < eps, &
                   msg="`idst(dst(x), 2)/(2.0_dp*(2.0_dp+1.0_dp))` failed.")
        call check(all(idst(dst(x, 2), 4)/(2.0_dp*(4.0_dp + 1.0_dp)) == &
                       [0.28138871112761987_dp, 0.78475214007354754_dp, 1.1919817084376492_dp, 0.94029999396468544_dp]), &
                   msg="`idst(dst(x, 2), 4)/(2.0_dp*(4.0_dp+1.0_dp))` failed.")

    end subroutine test_fftpack_idst

end program tester
