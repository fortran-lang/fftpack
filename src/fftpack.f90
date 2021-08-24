module fftpack

    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: zffti, zfftf, zfftb
    public :: fft, ifft
    public :: fftshift, ifftshift

    public :: dffti, dfftf, dfftb
    public :: rfft, irfft

    public :: dzffti, dzfftf, dzfftb

    public :: dcosqi, dcosqf, dcosqb
    public :: qct, iqct

    public :: dcosti, dcost
    public :: dct, idct

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
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine zfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `zfftf`.
        !> ([Specification](../page/specs/fftpack.html#zfftb))
        pure subroutine zfftb(n, c, wsave)
            import dp
            integer, intent(in) :: n
            complex(kind=dp), intent(inout) :: c(*)
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine zfftb

        !> Version: experimental
        !>
        !> Initialize `dfftf` and `dfftb`.
        !> ([Specification](../page/specs/fftpack.html#dffti))
        pure subroutine dffti(n, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(out) :: wsave(*)
        end subroutine dffti

        !> Version: experimental
        !>
        !> Forward transform of a double real periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#dfftf))
        pure subroutine dfftf(n, r, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(inout) :: r(*)
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine dfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dfftf`.
        !> ([Specification](../page/specs/fftpack.html#dfftb))
        pure subroutine dfftb(n, r, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(inout) :: r(*)
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine dfftb

        !> Version: experimental
        !>
        !> Initialize `dzfftf` and `dzfftb`.
        !> ([Specification](../page/specs/fftpack.html#dzffti))
        pure subroutine dzffti(n, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(out) :: wsave(*)
        end subroutine dzffti

        !> Version: experimental
        !>
        !> Simplified forward transform of a double real periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#dzfftf))
        pure subroutine dzfftf(n, r, azero, a, b, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(in) :: r(*)
            real(kind=dp), intent(out) :: azero
            real(kind=dp), intent(out) :: a(*), b(*)
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine dzfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dzfftf`.
        !> ([Specification](../page/specs/fftpack.html#dzfftb))
        pure subroutine dzfftb(n, r, azero, a, b, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(out) :: r(*)
            real(kind=dp), intent(in) :: azero
            real(kind=dp), intent(in) :: a(*), b(*)
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine dzfftb

        !> Version: experimental
        !>
        !> Initialize `dcosqf` and `dcosqb`.
        !> ([Specification](../page/specs/fftpack.html#dcosqi))
        pure subroutine dcosqi(n, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(out) :: wsave(*)
        end subroutine dcosqi

        !> Version: experimental
        !>
        !> Forward transform of quarter wave data.
        !> ([Specification](../page/specs/fftpack.html#dcosqf))
        pure subroutine dcosqf(n, x, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(inout) :: x(*)
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine dcosqf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dcosqf`.
        !> ([Specification](../page/specs/fftpack.html#dcosqb))
        pure subroutine dcosqb(n, x, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(inout) :: x(*)
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine dcosqb

        !> Version: experimental
        !>
        !> Initialize `dcost`. ([Specification](../page/specs/fftpack.html#dcosti))
        pure subroutine dcosti(n, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(out) :: wsave(*)
        end subroutine dcosti

        !> Version: experimental
        !>
        !> Discrete fourier cosine transform of an even sequence.
        !> ([Specification](../page/specs/fftpack.html#dcost))
        pure subroutine dcost(n, x, wsave)
            import dp
            integer, intent(in) :: n
            real(kind=dp), intent(inout) :: x(*)
            real(kind=dp), intent(in) :: wsave(*)
        end subroutine dcost

    end interface

    !> Version: experimental
    !>
    !> Forward transform of a double complex periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#fft))
    interface fft
        pure module function fft_dp(x, n) result(result)
            complex(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            complex(kind=dp), allocatable :: result(:)
        end function fft_dp
    end interface fft

    !> Version: experimental
    !>
    !> Backward transform of a double complex periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#ifft))
    interface ifft
        pure module function ifft_dp(x, n) result(result)
            complex(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            complex(kind=dp), allocatable :: result(:)
        end function ifft_dp
    end interface ifft

    !> Version: experimental
    !>
    !> Forward transform of a double real periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#rfft))
    interface rfft
        pure module function rfft_dp(x, n) result(result)
            real(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=dp), allocatable :: result(:)
        end function rfft_dp
    end interface rfft

    !> Version: experimental
    !>
    !> Backward transform of a double real periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#irfft))
    interface irfft
        pure module function irfft_dp(x, n) result(result)
            real(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=dp), allocatable :: result(:)
        end function irfft_dp
    end interface irfft

    !> Version: experimental
    !>
    !> Forward transform of quarter wave data.
    !> ([Specifiction](../page/specs/fftpack.html#qct))
    interface qct
        pure module function qct_dp(x, n) result(result)
            real(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=dp), allocatable :: result(:)
        end function qct_dp
    end interface qct

    !> Version: experimental
    !>
    !> Backward transform of quarter wave data.
    !> ([Specifiction](../page/specs/fftpack.html#iqct))
    interface iqct
        pure module function iqct_dp(x, n) result(result)
            real(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=dp), allocatable :: result(:)
        end function iqct_dp
    end interface iqct

    !> Version: experimental
    !>
    !> Discrete fourier cosine (forward) transform of an even sequence.
    !> ([Specification](../page/specs/fftpack.html#dct))
    interface dct
        pure module function dct_dp(x, n) result(result)
            real(kind=dp), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=dp), allocatable :: result(:)
        end function dct_dp
    end interface dct

    !> Version: experimental
    !>
    !> Discrete fourier cosine (backward) transform of an even sequence.
    !> ([Specification](../page/specs/fftpack.html#idct))
    interface idct
        module procedure :: dct_dp
    end interface idct

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
