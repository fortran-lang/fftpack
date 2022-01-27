module fftpack
    use fftpack_kind

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
    
    public :: rk

    interface

        !> Version: experimental
        !>
        !> Initialize `zfftf` and `zfftb`.
        !> ([Specification](../page/specs/fftpack.html#zffti))
        pure subroutine zffti(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine zffti

        !> Version: experimental
        !>
        !> Forward transform of a complex periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#zfftf))
        pure subroutine zfftf(n, c, wsave)
            import rk
            integer, intent(in) :: n
            complex(kind=rk), intent(inout) :: c(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine zfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `zfftf`.
        !> ([Specification](../page/specs/fftpack.html#zfftb))
        pure subroutine zfftb(n, c, wsave)
            import rk
            integer, intent(in) :: n
            complex(kind=rk), intent(inout) :: c(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine zfftb

        !> Version: experimental
        !>
        !> Initialize `dfftf` and `dfftb`.
        !> ([Specification](../page/specs/fftpack.html#dffti))
        pure subroutine dffti(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine dffti

        !> Version: experimental
        !>
        !> Forward transform of a real periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#dfftf))
        pure subroutine dfftf(n, r, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: r(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dfftf`.
        !> ([Specification](../page/specs/fftpack.html#dfftb))
        pure subroutine dfftb(n, r, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: r(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dfftb

        !> Version: experimental
        !>
        !> Initialize `dzfftf` and `dzfftb`.
        !> ([Specification](../page/specs/fftpack.html#dzffti))
        pure subroutine dzffti(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine dzffti

        !> Version: experimental
        !>
        !> Simplified forward transform of a real periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#dzfftf))
        pure subroutine dzfftf(n, r, azero, a, b, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(in) :: r(*)
            real(kind=rk), intent(out) :: azero
            real(kind=rk), intent(out) :: a(*), b(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dzfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dzfftf`.
        !> ([Specification](../page/specs/fftpack.html#dzfftb))
        pure subroutine dzfftb(n, r, azero, a, b, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: r(*)
            real(kind=rk), intent(in) :: azero
            real(kind=rk), intent(in) :: a(*), b(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dzfftb

        !> Version: experimental
        !>
        !> Initialize `dcosqf` and `dcosqb`.
        !> ([Specification](../page/specs/fftpack.html#dcosqi))
        pure subroutine dcosqi(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine dcosqi

        !> Version: experimental
        !>
        !> Forward transform of quarter wave data.
        !> ([Specification](../page/specs/fftpack.html#dcosqf))
        pure subroutine dcosqf(n, x, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: x(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dcosqf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dcosqf`.
        !> ([Specification](../page/specs/fftpack.html#dcosqb))
        pure subroutine dcosqb(n, x, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: x(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dcosqb

        !> Version: experimental
        !>
        !> Initialize `dcost`. ([Specification](../page/specs/fftpack.html#dcosti))
        pure subroutine dcosti(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine dcosti

        !> Version: experimental
        !>
        !> Discrete fourier cosine transform of an even sequence.
        !> ([Specification](../page/specs/fftpack.html#dcost))
        pure subroutine dcost(n, x, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: x(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dcost

    end interface

    !> Version: experimental
    !>
    !> Forward transform of a complex periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#fft))
    interface fft
        pure module function fft_rk(x, n) result(result)
            complex(kind=rk), intent(in) :: x(:)
            integer, intent(in), optional :: n
            complex(kind=rk), allocatable :: result(:)
        end function fft_rk
    end interface fft

    !> Version: experimental
    !>
    !> Backward transform of a complex periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#ifft))
    interface ifft
        pure module function ifft_rk(x, n) result(result)
            complex(kind=rk), intent(in) :: x(:)
            integer, intent(in), optional :: n
            complex(kind=rk), allocatable :: result(:)
        end function ifft_rk
    end interface ifft

    !> Version: experimental
    !>
    !> Forward transform of a real periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#rfft))
    interface rfft
        pure module function rfft_rk(x, n) result(result)
            real(kind=rk), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=rk), allocatable :: result(:)
        end function rfft_rk
    end interface rfft

    !> Version: experimental
    !>
    !> Backward transform of a real periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#irfft))
    interface irfft
        pure module function irfft_rk(x, n) result(result)
            real(kind=rk), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=rk), allocatable :: result(:)
        end function irfft_rk
    end interface irfft

    !> Version: experimental
    !>
    !> Forward transform of quarter wave data.
    !> ([Specifiction](../page/specs/fftpack.html#qct))
    interface qct
        pure module function qct_rk(x, n) result(result)
            real(kind=rk), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=rk), allocatable :: result(:)
        end function qct_rk
    end interface qct

    !> Version: experimental
    !>
    !> Backward transform of quarter wave data.
    !> ([Specifiction](../page/specs/fftpack.html#iqct))
    interface iqct
        pure module function iqct_rk(x, n) result(result)
            real(kind=rk), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=rk), allocatable :: result(:)
        end function iqct_rk
    end interface iqct

    !> Version: experimental
    !>
    !> Discrete fourier cosine (forward) transform of an even sequence.
    !> ([Specification](../page/specs/fftpack.html#dct))
    interface dct
        pure module function dct_rk(x, n) result(result)
            real(kind=rk), intent(in) :: x(:)
            integer, intent(in), optional :: n
            real(kind=rk), allocatable :: result(:)
        end function dct_rk
    end interface dct

    !> Version: experimental
    !>
    !> Discrete fourier cosine (backward) transform of an even sequence.
    !> ([Specification](../page/specs/fftpack.html#idct))
    interface idct
        module procedure :: dct_rk
    end interface idct

    !> Version: experimental
    !>
    !> Shifts zero-frequency component to center of spectrum.
    !> ([Specifiction](../page/specs/fftpack.html#fftshift))
    interface fftshift
        pure module function fftshift_crk(x) result(result)
            complex(kind=rk), intent(in) :: x(:)
            complex(kind=rk), allocatable :: result(:)
        end function fftshift_crk
        pure module function fftshift_rrk(x) result(result)
            real(kind=rk), intent(in) :: x(:)
            real(kind=rk), allocatable :: result(:)
        end function fftshift_rrk
    end interface fftshift

    !> Version: experimental
    !>
    !> Shifts zero-frequency component to beginning of spectrum.
    !> ([Specifiction](../page/specs/fftpack.html#ifftshift))
    interface ifftshift
        pure module function ifftshift_crk(x) result(result)
            complex(kind=rk), intent(in) :: x(:)
            complex(kind=rk), allocatable :: result(:)
        end function ifftshift_crk
        pure module function ifftshift_rrk(x) result(result)
            real(kind=rk), intent(in) :: x(:)
            real(kind=rk), allocatable :: result(:)
        end function ifftshift_rrk
    end interface ifftshift

end module fftpack
