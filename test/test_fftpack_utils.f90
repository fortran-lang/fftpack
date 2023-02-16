module test_fftpack_utils

    use fftpack, only: rk, fft, ifft, fftshift, ifftshift, fftfreq, rfftfreq
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
                    new_unittest("ifftshift_real", test_fftshift_real), &
                    new_unittest("fftfreq_1", test_fftfreq_1), &
                    new_unittest("fftfreq_2", test_fftfreq_2), &
                    new_unittest("fftfreq_3", test_fftfreq_3), &
                    new_unittest("rfftfreq", test_rfftfreq) &
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

    subroutine test_fftfreq_1(error)
        type(error_type), allocatable, intent(out) :: error
        integer, dimension(8) :: target1 = [0, 1, 2, 3, -4, -3, -2, -1]
        integer, dimension(9) :: target2 = [0, 1, 2, 3, 4, -4, -3, -2, -1]

        call check(error, all(fftfreq(8) == target1),&
                    "all(fftfreq(8) == target1) failed.")
        if (allocated(error)) return
        call check(error, all(fftfreq(9) == target2),&
                    "all(fftfreq(9) == target2) failed.")
    end subroutine test_fftfreq_1

    subroutine test_fftfreq_2(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error

        real(rk), parameter :: tol = 1.0e-12_rk
        real(rk), parameter :: twopi = 8*atan(1.0_rk)  ! 2*pi
        complex(rk), parameter :: imu = (0,1)  ! imaginary unit

        integer, parameter :: n = 128
        integer :: i
        complex(rk), dimension(n) :: xvec, xfou
        real(rk), dimension(n) :: xtrue

        do  i = 1, n
            xvec(i) = cos(twopi*(i-1)/n)
            xtrue(i) = -sin(twopi*(i-1)/n)  ! derivative in physical space
        end do

        xfou = fft(xvec)/n
        xfou = imu*fftfreq(n)*xfou  ! derivative in Fourier space
        xvec = ifft(xfou)
        call check(error, maxval(abs(xvec-xtrue)) < tol, &
                    "maxval(abs(xvec-xtrue)) < tol failed.")
    end subroutine test_fftfreq_2

    subroutine test_fftfreq_3(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error

        real(rk), parameter :: tol = 1.0e-12_rk
        real(rk), parameter :: twopi = 8*atan(1.0_rk)  ! 2*pi
        complex(rk), parameter :: imu = (0,1)  ! imaginary unit

        integer, parameter :: n = 135
        integer :: i
        complex(rk), dimension(n) :: xvec, xfou
        real(rk), dimension(n) :: xtrue

        do  i = 1, n
            xvec(i) = cos(twopi*(i-1)/n)
            xtrue(i) = -sin(twopi*(i-1)/n)  ! derivative in physical space
        end do

        xfou = fft(xvec)/n
        xfou = imu*fftfreq(n)*xfou  ! derivative in Fourier space
        xvec = ifft(xfou)
        call check(error, maxval(abs(xvec-xtrue)) < tol, &
                    "maxval(abs(xvec-xtrue)) < tol failed.")
    end subroutine test_fftfreq_3

    subroutine test_rfftfreq(error)
        type(error_type), allocatable, intent(out) :: error
        integer, dimension(8) :: target1 = [0, 1, 1, 2, 2, 3, 3, -4]
        integer, dimension(9) :: target2 = [0, 1, 1, 2, 2, 3, 3, 4, 4]

        call check(error, all(rfftfreq(8) == target1),&
                    "all(rfftfreq(8) == target1) failed.")
        if (allocated(error)) return
        call check(error, all(rfftfreq(9) == target2),&
                    "all(rfftfreq(9) == target2) failed.")
    end subroutine test_rfftfreq

end module test_fftpack_utils
