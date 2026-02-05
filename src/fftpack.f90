module fftpack
   use fftpack_kind, only: rk

   implicit none
   private

   public :: zffti, zfftf, zfftb
   public :: fft, ifft
   public :: fftshift, ifftshift
   public :: fftfreq, rfftfreq

   public :: dffti, dfftf, dfftb
   public :: rfft, irfft

   public :: dzffti, dzfftf, dzfftb

   public :: dcosqi, dcosqf, dcosqb
   public :: dcosti, dcost
   public :: dct, idct
   public :: dct_t1i, dct_t1
   public :: dct_t23i, dct_t2, dct_t3

   public :: dsinti, dsint

   public :: rk

   interface

      !> Version: experimental
      !>
      !> Initialize `zfftf` and `zfftb`.
      !> ([Specification](../page/specs/fftpack.html#zffti))
      pure subroutine zffti(n, wsave)
         import rk
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(out) :: wsave(*)
      end subroutine zffti

      !> Version: experimental
      !>
      !> Forward transform of a complex periodic sequence.
      !> ([Specification](../page/specs/fftpack.html#zfftf))
      pure subroutine zfftf(n, c, wsave)
         import rk
         implicit none
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
         implicit none
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
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(out) :: wsave(*)
      end subroutine dffti

      !> Version: experimental
      !>
      !> Forward transform of a real periodic sequence.
      !> ([Specification](../page/specs/fftpack.html#dfftf))
      pure subroutine dfftf(n, r, wsave)
         import rk
         implicit none
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
         implicit none
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
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(out) :: wsave(*)
      end subroutine dzffti

      !> Version: experimental
      !>
      !> Simplified forward transform of a real periodic sequence.
      !> ([Specification](../page/specs/fftpack.html#dzfftf))
      pure subroutine dzfftf(n, r, azero, a, b, wsave)
         import rk
         implicit none
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
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(out) :: r(*)
         real(kind=rk), intent(in) :: azero
         real(kind=rk), intent(in) :: a(*), b(*)
         real(kind=rk), intent(in) :: wsave(*)
      end subroutine dzfftb

      !> Version: experimental
      !>
      !> Initialize `dcosqf` and `dcosqb`.
      !> ([Specification](../page/specs/fftpack.html#initialize-dct-2-3-dcosqi-or-dct_t23i))
      pure subroutine dcosqi(n, wsave)
         import rk
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(out) :: wsave(*)
      end subroutine dcosqi

      !> Version: experimental
      !>
      !> Forward transform of quarter wave data.
      !> ([Specification](../page/specs/fftpack.html#compute-dct-3-dcosqf-or-dct_t3))
      pure subroutine dcosqf(n, x, wsave)
         import rk
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(inout) :: x(*)
         real(kind=rk), intent(in) :: wsave(*)
      end subroutine dcosqf

      !> Version: experimental
      !>
      !> Unnormalized inverse of `dcosqf`.
      !> ([Specification](../page/specs/fftpack.html#compute-dct-2-dcosqb-or-dct_t2))
      pure subroutine dcosqb(n, x, wsave)
         import rk
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(inout) :: x(*)
         real(kind=rk), intent(in) :: wsave(*)
      end subroutine dcosqb

      !> Version: experimental
      !>
      !> Initialize `dcost`.
      !> ([Specification](../page/specs/fftpack.html#initialize-dct-1-dcosti-or-dct_t1i))
      pure subroutine dcosti(n, wsave)
         import rk
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(out) :: wsave(*)
      end subroutine dcosti

      !> Version: experimental
      !>
      !> Discrete fourier cosine transform of an even sequence.
      !> ([Specification](../page/specs/fftpack.html#compute-dct-1-dcost-or-dct_t1))
      pure subroutine dcost(n, x, wsave)
         import rk
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(inout) :: x(*)
         real(kind=rk), intent(in) :: wsave(*)
      end subroutine dcost

      pure subroutine dsinti(n, wsave)
         import rk
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(out) :: wsave(*)
      end subroutine dsinti

      pure subroutine dsint(n, x, wsave)
         import rk
         implicit none
         integer, intent(in) :: n
         real(kind=rk), intent(inout) :: x(*)
         real(kind=rk), intent(in) :: wsave(*)
      end subroutine dsint

      !> Version: experimental
      !>
      !> Integer frequency values involved in complex FFT.
      !> ([Specifiction](../page/specs/fftpack.html#fftfreq))
      pure module function fftfreq(n) result(out)
         implicit none
         integer, intent(in) :: n
         integer, dimension(n) :: out
      end function fftfreq

      !> Version: experimental
      !>
      !> Integer frequency values involved in real FFT.
      !> ([Specifiction](../page/specs/fftpack.html#rfftfreq))
      pure module function rfftfreq(n) result(out)
         implicit none
         integer, intent(in) :: n
         integer, dimension(n) :: out
      end function rfftfreq

   end interface

   !> Version: experimental
   !>
   !> Forward transform of a complex periodic sequence.
   !> ([Specifiction](../page/specs/fftpack.html#fft))
   interface fft
      pure module function fft_rk(x, n) result(result)
         implicit none
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
         implicit none
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
         implicit none
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
         implicit none
         real(kind=rk), intent(in) :: x(:)
         integer, intent(in), optional :: n
         real(kind=rk), allocatable :: result(:)
      end function irfft_rk
   end interface irfft

   !> Version: experimental
   !>
   !> Dsicrete cosine transforms.
   !> ([Specification](../page/specs/fftpack.html#simplified-dct-of-types-1-2-3-dct))
   interface dct
      pure module function dct_rk(x, n, type) result(result)
         implicit none
         real(kind=rk), intent(in) :: x(:)
         integer, intent(in), optional :: n
         integer, intent(in), optional :: type
         real(kind=rk), allocatable :: result(:)
      end function dct_rk
   end interface dct

   !> Version: experimental
   !>
   !> Inverse discrete cosine transforms.
   !> ([Specification](../page/specs/fftpack.html#simplified-inverse-dct-of-types-1-2-3-idct))
   interface idct
      pure module function idct_rk(x, n, type) result(result)
         implicit none
         real(kind=rk), intent(in) :: x(:)
         integer, intent(in), optional :: n
         integer, intent(in), optional :: type
         real(kind=rk), allocatable :: result(:)
      end function idct_rk
   end interface idct

   !> Version: experimental
   !>
   !> Initialize DCT type-1
   !> ([Specification](../page/specs/fftpack.html#initialize-dct-1-dcosti-or-dct_t1i))
   interface dct_t1i
      procedure :: dcosti
   end interface dct_t1i

   !> Version: experimental
   !>
   !> Perform DCT type-1
   !> ([Specification](../page/specs/fftpack.html#compute-dct-1-dcost-or-dct_t1))
   interface dct_t1
      procedure :: dcost
   end interface dct_t1

   !> Version: experimental
   !>
   !> Initialize DCT types 2, 3
   !> ([Specification](../page/specs/fftpack.html#initialize-dct-2-3-dcosqi-or-dct_t23i))
   interface dct_t23i
      procedure :: dcosqi
   end interface dct_t23i

   !> Version: experimental
   !>
   !> Perform DCT type-2
   !> ([Specification](../page/specs/fftpack.html#compute-dct-2-dcosqb-or-dct_t2))
   interface dct_t2
      procedure :: dcosqb
   end interface dct_t2

   !> Version: experimental
   !>
   !> Perform DCT type-3
   !> ([Specification](../page/specs/fftpack.html#compute-dct-3-dcosqf-or-dct_t3))
   interface dct_t3
      procedure :: dcosqf
   end interface dct_t3

   !> Version: experimental
   !>
   !> Shifts zero-frequency component to center of spectrum.
   !> ([Specifiction](../page/specs/fftpack.html#fftshift))
   interface fftshift
      pure module function fftshift_crk(x) result(result)
         implicit none
         complex(kind=rk), intent(in) :: x(:)
         complex(kind=rk), dimension(size(x)) :: result
      end function fftshift_crk
      pure module function fftshift_rrk(x) result(result)
         implicit none
         real(kind=rk), intent(in) :: x(:)
         real(kind=rk), dimension(size(x)) :: result
      end function fftshift_rrk
   end interface fftshift

   !> Version: experimental
   !>
   !> Shifts zero-frequency component to beginning of spectrum.
   !> ([Specifiction](../page/specs/fftpack.html#ifftshift))
   interface ifftshift
      pure module function ifftshift_crk(x) result(result)
         implicit none
         complex(kind=rk), intent(in) :: x(:)
         complex(kind=rk), dimension(size(x)) :: result
      end function ifftshift_crk
      pure module function ifftshift_rrk(x) result(result)
         implicit none
         real(kind=rk), intent(in) :: x(:)
         real(kind=rk), dimension(size(x)) :: result
      end function ifftshift_rrk
   end interface ifftshift

end module fftpack
