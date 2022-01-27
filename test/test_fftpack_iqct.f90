program tester

    call test_fftpack_iqct()
    print *, "All tests in `test_fftpack_iqct` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_iqct
        use fftpack, only: qct, iqct
    use fftpack_kind
        real(kind=rk) :: eps = 1.0e-10_rk

        real(kind=rk) :: x(4) = [1, 2, 3, 4]

        call check(sum(abs(iqct(qct(x))/(4.0_rk*4.0_rk) - [real(kind=rk) :: 1, 2, 3, 4])) < eps, &
                   msg="`iqct(qct(x)/(4.0_rk*4.0_rk)` failed.")
        call check(sum(abs(iqct(qct(x), 2)/(4.0_rk*2.0_rk) - [1.4483415291679655_rk, 7.4608849947753271_rk])) < eps, &
                   msg="`iqct(qct(x), 2)/(4.0_rk*2.0_rk)` failed.")
        call check(sum(abs(iqct(qct(x, 2), 4)/(4.0_rk*4.0_rk) - [0.5_rk, 0.70932417358418376_rk, 1.0_rk, &
                                                                 0.78858050747473762_rk])) < eps, &
                   msg="`iqct(qct(x, 2), 4)/(4.0_rk*4.0_rk)` failed.")

    end subroutine test_fftpack_iqct

end program tester
