module test_fftpack_utils

    use fftpack, only: rk, fftshift, ifftshift
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_utils

contains

    !> Collect all exported unit tests
    subroutine collect_utils(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("fftshift_complex", test_fftshift_complex), &
                    new_unittest("fftshift_real", test_fftshift_real), &
                    new_unittest("ifftshift_complex", test_fftshift_complex), &
                    new_unittest("ifftshift_real", test_fftshift_real) &
                    ]

    end subroutine collect_utils

    subroutine test_fftshift_complex(error)
        type(error_type), allocatable, intent(out) :: error
        complex(kind=rk) :: xeven(4) = [1, 2, 3, 4]
        complex(kind=rk) :: xodd(5) = [1, 2, 3, 4, 5]

        call check(error, all(fftshift(xeven) == [complex(kind=rk) :: 3, 4, 1, 2]), &
                   "all(fftshift(xeven) == [complex(kind=rk) :: 3, 4, 1, 2]) failed.")
        if (allocated(error)) return
        call check(error, all(fftshift(xodd) == [complex(kind=rk) :: 4, 5, 1, 2, 3]), &
                   "all(fftshift(xodd) == [complex(kind=rk) :: 4, 5, 1, 2, 3]) failed.")

    end subroutine test_fftshift_complex

    subroutine test_fftshift_real(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: xeven(4) = [1, 2, 3, 4]
        real(kind=rk) :: xodd(5) = [1, 2, 3, 4, 5]

        call check(error, all(fftshift(xeven) == [real(kind=rk) :: 3, 4, 1, 2]), &
                   "all(fftshift(xeven) == [real(kind=rk) :: 3, 4, 1, 2]) failed.")
        if (allocated(error)) return
        call check(error, all(fftshift(xodd) == [real(kind=rk) :: 4, 5, 1, 2, 3]), &
                   "all(fftshift(xodd) == [real(kind=rk) :: 4, 5, 1, 2, 3]) failed.")

    end subroutine test_fftshift_real

    subroutine test_ifftshift_complex(error)
        type(error_type), allocatable, intent(out) :: error
        integer :: i

        complex(kind=rk) :: xeven(4) = [3, 4, 1, 2]
        complex(kind=rk) :: xodd(5) = [4, 5, 1, 2, 3]

        call check(error, all(ifftshift(xeven) == [complex(kind=rk) ::(i, i=1, 4)]), &
                   "all(ifftshift(xeven) == [complex(kind=rk) ::(i, i=1, 4)]) failed.")
        if (allocated(error)) return
        call check(error, all(ifftshift(xodd) == [complex(kind=rk) ::(i, i=1, 5)]), &
                   "all(ifftshift(xodd) == [complex(kind=rk) ::(i, i=1, 5)]) failed.")

    end subroutine test_ifftshift_complex

    subroutine test_ifftshift_real(error)
        type(error_type), allocatable, intent(out) :: error
        integer :: i

        real(kind=rk) :: xeven(4) = [3, 4, 1, 2]
        real(kind=rk) :: xodd(5) = [4, 5, 1, 2, 3]

        call check(error, all(ifftshift(xeven) == [real(kind=rk) ::(i, i=1, 4)]), &
                   "all(ifftshift(xeven) == [real(kind=rk) ::(i, i=1, 4)]) failed.")
        if (allocated(error)) return
        call check(error, all(ifftshift(xodd) == [real(kind=rk) ::(i, i=1, 5)]), &
                   "all(ifftshift(xodd) == [real(kind=rk) ::(i, i=1, 5)]) failed.")

    end subroutine test_ifftshift_real

end module test_fftpack_utils
