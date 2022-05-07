module test_fftpack_qct

    use fftpack, only: rk, dcosqi, dcosqf, dcosqb, qct, iqct
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_qct

#if defined(fftpack_sp)
    real(kind=rk) :: eps = 1.0e-5_rk
#else
    real(kind=rk) :: eps = 1.0e-10_rk
#endif

contains

    !> Collect all exported unit tests
    subroutine collect_qct(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("classic-qct-API", test_classic_qct), &
                    new_unittest("modernized-qct-API", test_modernized_qct), &
                    new_unittest("modernized-iqct-API", test_modernized_iqct) &
                    ]

    end subroutine collect_qct

    subroutine test_classic_qct(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: w(3*4 + 15)
        real(kind=rk) :: x(4) = [1, 2, 3, 4]

        call dcosqi(4, w)
        call dcosqf(4, x, w)
        call check(error, sum(abs(x - [11.999626276085150_rk, -9.1029432177492193_rk, &
                                       2.6176618435106480_rk, -1.5143449018465791_rk])) < eps, &
                   "`dcosqf` failed.")
        if (allocated(error)) return
        call dcosqb(4, x, w)
        call check(error, sum(abs(x/(4.0_rk*4.0_rk) - [real(kind=rk) :: 1, 2, 3, 4])) < eps, &
                   "`dcosqb` failed.")

    end subroutine test_classic_qct

    subroutine test_modernized_qct(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: x(3) = [9, -9, 3]

        call check(error, sum(abs(qct(x, 2) - [-3.7279220613578570_rk, 21.727922061357859_rk])) < eps, &
                   "`qct(x, 2)` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(qct(x, 3) - qct(x))) < eps, &
                   "`qct(x,3)` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(qct(x, 4) - [-3.3871908980838743_rk, -2.1309424696909023_rk, &
                                               11.645661095452331_rk, 29.872472272322447_rk])) < eps, &
                   "`qct(x, 4)` failed.")

    end subroutine test_modernized_qct

    subroutine test_modernized_iqct(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: x(4) = [1, 2, 3, 4]

        call check(error, sum(abs(iqct(qct(x))/(4.0_rk*4.0_rk) - [real(kind=rk) :: 1, 2, 3, 4])) < eps, &
                   "`iqct(qct(x)/(4.0_rk*4.0_rk)` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(iqct(qct(x), 2)/(4.0_rk*2.0_rk) - [1.4483415291679655_rk, 7.4608849947753271_rk])) < eps, &
                   "`iqct(qct(x), 2)/(4.0_rk*2.0_rk)` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(iqct(qct(x, 2), 4)/(4.0_rk*4.0_rk) - [0.5_rk, 0.70932417358418376_rk, 1.0_rk, &
                                                                        0.78858050747473762_rk])) < eps, &
                   "`iqct(qct(x, 2), 4)/(4.0_rk*4.0_rk)` failed.")

    end subroutine test_modernized_iqct

end module test_fftpack_qct
