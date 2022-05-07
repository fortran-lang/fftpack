module test_fftpack_fft

    use fftpack, only: rk, zffti, zfftf, zfftb, fft, ifft
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_fft

contains

    !> Collect all exported unit tests
    subroutine collect_fft(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("classic-fft-API", test_classic_fft), &
                    new_unittest("modernized-fft-API", test_modernized_fft), &
                    new_unittest("modernized-ifft-API", test_modernized_ifft) &
                    ]

    end subroutine collect_fft

    subroutine test_classic_fft(error)
        type(error_type), allocatable, intent(out) :: error
        complex(kind=rk) :: x(4) = [1, 2, 3, 4]
        real(kind=rk) :: w(31)

        call zffti(4, w)
        call zfftf(4, x, w)
        call check(error, all(x == [complex(kind=rk) ::(10, 0), (-2, 2), (-2, 0), (-2, -2)]), &
                   "`zfftf` failed.")
        if (allocated(error)) return
        call zfftb(4, x, w)
        call check(error, all(x/4.0_rk == [complex(kind=rk) ::(1, 0), (2, 0), (3, 0), (4, 0)]), &
                   "`zfftb` failed.")

    end subroutine test_classic_fft

    subroutine test_modernized_fft(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: eps = 1.0e-10_rk
        complex(kind=rk) :: x(3) = [1.0_rk, 2.0_rk, 3.0_rk]

        call check(error, sum(abs(fft(x, 2) - [(3.0_rk, 0.0_rk), (-1.0_rk, 0.0_rk)])) < eps, &
                   "`fft(x, 2)` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(fft(x, 3) - fft(x))) < eps, &
                   "`fft(x, 3)` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(fft(x, 4) - [(6.0_rk, 0.0_rk), (-2.0_rk, -2.0_rk), (2.0_rk, 0.0_rk), (-2.0_rk, 2.0_rk)])) < eps, &
                   "`fft(x, 4)` failed.")

    end subroutine test_modernized_fft
    
    subroutine test_modernized_ifft(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=rk) :: eps = 1.0e-10_rk
        complex(kind=rk) :: x(4) = [1, 2, 3, 4]
        
        call check(error, sum(abs(ifft(fft(x))/4.0_rk - [complex(kind=rk) :: 1, 2, 3, 4])) < eps, &
                   "`ifft(fft(x))/4.0_rk` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(ifft(fft(x), 2) - [complex(kind=rk) ::(8, 2), (12, -2)])) < eps, &
                   "`ifft(fft(x), 2)` failed.")
        if (allocated(error)) return
        call check(error, sum(abs(ifft(fft(x, 2), 4) - [complex(kind=rk) ::(2, 0), (3, -1), (4, 0), (3, 1)])) < eps, &
                   "`ifft(fft(x, 2), 4)` failed.")

    end subroutine test_modernized_ifft
    

end module test_fftpack_fft
