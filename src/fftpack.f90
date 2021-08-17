module fftpack

    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: zffti, zfftf, zfftb
    public :: fft, ifft
    public :: fftshift, ifftshift

    interface

        !> Version: experimental
        !>
        !> Initialize `zfftf` and `zfftb`.
        !> ([Specification](../page/specs/fftpack.html#zffti))
        pure subroutine zffti(n, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(out) :: wsave(*)
        end subroutine zffti

        !> Version: experimental
        !>
        !> Forward transform of a double complex periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#zfftf))
        pure subroutine zfftf(n, c, wsave)
            import dp
            integer, intent(in) :: n
            complex(kind=dp), intent(inout) :: c(*)
            real(kind=dp), intent(inout) :: wsave(*)
        end subroutine zfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `zfftf`.
        !> ([Specification](../page/specs/fftpack.html#zfftb))
        pure subroutine zfftb(n, c, wsave)
            import dp
            integer, intent(in) :: n
            complex(kind=dp), intent(inout) :: c(*)
            real(kind=dp), intent(inout) :: wsave(*)
        end subroutine zfftb

    end interface

    !> Version: experimental
    !>
    !> Forward transform of a double complex periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#fft))
    interface fft
        pure module function fft_cdp(x, n) result(result)
            complex(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            complex(kind=dp), allocatable :: result(:)
        end function fft_cdp
    end interface fft

    !> Version: experimental
    !>
    !> Backward transform of a double complex periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#ifft))
    interface ifft
        pure module function ifft_cdp(x, n) result(result)
            complex(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            complex(kind=dp), allocatable :: result(:)
        end function ifft_cdp
    end interface ifft

    !> Version: experimental
    !>
    !> Shifts zero-frequency component to center of spectrum.
    !> ([Specifiction](../page/specs/fftpack.html#fftshift))
    interface fftshift
        pure module function fftshift_cdp(x) result(result)
            complex(kind=dp), intent(in) :: x(:)
            complex(kind=dp), allocatable :: result(:)
        end function fftshift_cdp
        pure module function fftshift_rdp(x) result(result)
            real(kind=dp), intent(in) :: x(:)
            real(kind=dp), allocatable :: result(:)
        end function fftshift_rdp
    end interface fftshift

    !> Version: experimental
    !>
    !> Shifts zero-frequency component to beginning of spectrum.
    !> ([Specifiction](../page/specs/fftpack.html#ifftshift))
    interface ifftshift
        pure module function ifftshift_cdp(x) result(result)
            complex(kind=dp), intent(in) :: x(:)
            complex(kind=dp), allocatable :: result(:)
        end function ifftshift_cdp
        pure module function ifftshift_rdp(x) result(result)
            real(kind=dp), intent(in) :: x(:)
            real(kind=dp), allocatable :: result(:)
        end function ifftshift_rdp
    end interface ifftshift

end module fftpack
