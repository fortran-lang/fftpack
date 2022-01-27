program tester

    call test_fftpack_dcosq_real
    print *, "All tests in `test_fftpack_dcosq` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dcosq_real
        use fftpack, only: dcosqi, dcosqf, dcosqb
        use fftpack_kind
        real(kind=rk) :: w(3*4 + 15)
        real(kind=rk) :: x(4) = [1, 2, 3, 4]
        real(kind=rk) :: eps = 1.0e-10_rk

        call dcosqi(4, w)
        call dcosqf(4, x, w)
        call check(sum(abs(x - [11.999626276085150_rk, -9.1029432177492193_rk, &
                                2.6176618435106480_rk, -1.5143449018465791_rk])) < eps, msg="`dcosqf` failed.")
        call dcosqb(4, x, w)
        call check(sum(abs(x/(4.0_rk*4.0_rk) - [real(kind=rk) :: 1, 2, 3, 4])) < eps, msg="`dcosqb` failed.")

    end subroutine test_fftpack_dcosq_real

end program tester
