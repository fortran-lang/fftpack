module test_fftpack_dct

    use fftpack, only: rk, dcosti, dcost, dct, idct
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_dct

#if defined(fftpack_sp)
    real(kind=rk) :: eps = 1.0e-5_rk
#else
    real(kind=rk) :: eps = 1.0e-10_rk
#endif

contains

    !> Collect all exported unit tests
    subroutine collect_dct(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("classic-dct-API", test_classic_dct), &
                    new_unittest("modernized-dct-API", test_modernized_dct), &
                    new_unittest("modernized-idct-API", test_modernized_idct) &
                    ]

    end subroutine collect_dct

    subroutine test_classic_dct(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: w(3*4 + 15)
        real(kind=rk) :: x(4) = [1, 2, 3, 4]

        call dcosti(4, w)
        call dcost(4, x, w)
        call check(error, sum(abs(x - [real(kind=rk) :: 15, -4, 0, -1.0000000000000009_rk])) < eps, &
                   "`dcosti` failed.")
        if (allocated(error)) return

        call dcost(4, x, w)
        call check(error, sum(abs(x/(2.0_rk*(4.0_rk - 1.0_rk)) - &
                                  [real(kind=rk) :: 1, 2, 3, 4])) < eps, "`dcost` failed.")

    end subroutine test_classic_dct

    subroutine test_modernized_dct(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: x(3) = [9, -9, 3]

        call check(error, all(dct(x, 2) == [real(kind=rk) :: 0, 18]), "`dct(x, 2)` failed.")
        if (allocated(error)) return
        call check(error, all(dct(x, 3) == dct(x)), "`dct(x, 3)` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(dct(x, 4) - [real(kind=rk) :: -3, -3.0000000000000036_rk, 15, 33])) &
                   < eps, "`dct(x, 4)` failed.")

    end subroutine test_modernized_dct

    subroutine test_modernized_idct(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: x(4) = [1, 2, 3, 4]

        call check(error, sum(abs(idct(dct(x))/(2.0_rk*(4.0_rk - 1.0_rk)) - &
                                  [real(kind=rk) :: 1, 2, 3, 4])) < eps, &
                   "`idct(dct(x))/(2.0_rk*(4.0_rk-1.0_rk))` failed.")
        if (allocated(error)) return
        call check(error, all(idct(dct(x), 2)/(2.0_rk*(2.0_rk - 1.0_rk)) == [real(kind=rk) :: 5.5, 9.5]), &
                   "`idct(dct(x), 2)/(2.0_rk*(2.0_rk-1.0_rk))` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(idct(dct(x, 2), 4)/(2.0_rk*(4.0_rk - 1.0_rk)) - &
                                  [0.16666666666666666_rk, 0.33333333333333331_rk, &
                                   0.66666666666666663_rk, 0.83333333333333315_rk])) < eps, &
                   "`idct(dct(x, 2), 4)/(2.0_rk*(4.0_rk-1.0_rk))` failed.")

    end subroutine test_modernized_idct

end module test_fftpack_dct
