submodule(fftpack) fftpack_fftshift
   implicit none
contains

   !> Shifts zero-frequency component to center of spectrum for `complex` type.
   pure module function fftshift_cdp(x) result(result)
      complex(kind=dp), intent(in) :: x(:)
      complex(kind=dp), dimension(size(x)) :: result

      result = cshift(x, shift=-floor(0.5_dp*size(x)))

   end function fftshift_cdp

   !> Shifts zero-frequency component to center of spectrum for `real` type.
   pure module function fftshift_rdp(x) result(result)
      real(kind=dp), intent(in) :: x(:)
      real(kind=dp), dimension(size(x)) :: result

      result = cshift(x, shift=-floor(0.5_dp*size(x)))

   end function fftshift_rdp

end submodule fftpack_fftshift
