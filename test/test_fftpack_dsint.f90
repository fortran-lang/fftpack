program tester

    call test_fftpack_dsint_real
    print *, "All tests in `test_fftpack_dsint` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dsint_real
        use fftpack, only: dsinti, dsint
        use iso_fortran_env, only: dp => real64
        real(kind=dp) :: w(3*4 + 15)
        real(kind=dp) :: x(4) = [1, 2, 3, 4]
        real(kind=dp) :: eps = 1.0e-10_dp

        call dsinti(4, w)
        call dsint(4, x, w)
        call check(sum(abs(x - [15.388417685876266_dp, -6.8819096023558668_dp, &
                                3.6327126400268046_dp, -1.6245984811645318_dp])) < eps, &
                   msg="`dsinti` failed.")

        call dsint(4, x, w)
        call check(sum(abs(x/(2.0_dp*(4.0_dp + 1.0_dp)) - [real(kind=dp) :: 1, 2, 3, 4])) < eps, &
                   msg="`dsint` failed.")

    end subroutine test_fftpack_dsint_real

end program tester
