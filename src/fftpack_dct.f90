submodule(fftpack) fftpack_dct
   implicit none
contains

   !> Discrete cosine transforms of types 1, 2, 3.
   pure module function dct_dp(x, n, type) result(result)
      real(kind=dp), intent(in) :: x(:)
      integer, intent(in), optional :: n
      integer, intent(in), optional :: type
      real(kind=dp), allocatable :: result(:)

      integer :: lenseq, lensav, i
      real(kind=dp), allocatable :: wsave(:)

      if (present(n)) then
         lenseq = n
         if (lenseq <= size(x)) then
            result = x(:lenseq)
         else if (lenseq > size(x)) then
            result = [x, (0.0_dp, i=1, lenseq - size(x))]
         end if
      else
         lenseq = size(x)
         result = x
      end if

      ! Default to DCT-2
      if (.not. present(type)) then
         lensav = 3*lenseq + 15
         allocate (wsave(lensav))
         call dcosqi(lenseq, wsave)
         call dcosqb(lenseq, result, wsave)
         return
      end if

      if (type == 1) then  ! DCT-1
         lensav = 3*lenseq + 15
         allocate (wsave(lensav))
         call dcosti(lenseq, wsave)
         call dcost(lenseq, result, wsave)
      else if (type == 2) then  ! DCT-2
         lensav = 3*lenseq + 15
         allocate (wsave(lensav))
         call dcosqi(lenseq, wsave)
         call dcosqb(lenseq, result, wsave)
      else if (type == 3) then  ! DCT-3
         lensav = 3*lenseq + 15
         allocate (wsave(lensav))
         call dcosqi(lenseq, wsave)
         call dcosqf(lenseq, result, wsave)
      end if
   end function dct_dp

   !> Inverse discrete cosine transforms of types 1, 2, 3.
   pure module function idct_dp(x, n, type) result(result)
      real(kind=dp), intent(in) :: x(:)
      integer, intent(in), optional :: n
      integer, intent(in), optional :: type
      real(kind=dp), allocatable :: result(:)

      integer :: lenseq, lensav, i
      real(kind=dp), allocatable :: wsave(:)

      if (present(n)) then
         lenseq = n
         if (lenseq <= size(x)) then
            result = x(:lenseq)
         else if (lenseq > size(x)) then
            result = [x, (0.0_dp, i=1, lenseq - size(x))]
         end if
      else
         lenseq = size(x)
         result = x
      end if

      ! Default to t=2; inverse DCT-2 is DCT-3
      if (.not. present(type)) then
         lensav = 3*lenseq + 15
         allocate (wsave(lensav))
         call dcosqi(lenseq, wsave)
         call dcosqf(lenseq, result, wsave)
         return
      end if

      if (type == 1) then  ! inverse DCT-1 is DCT-1
         lensav = 3*lenseq + 15
         allocate (wsave(lensav))
         call dcosti(lenseq, wsave)
         call dcost(lenseq, result, wsave)
      else if (type == 2) then  ! inverse DCT-2 is DCT-3
         lensav = 3*lenseq + 15
         allocate (wsave(lensav))
         call dcosqi(lenseq, wsave)
         call dcosqf(lenseq, result, wsave)
      else if (type == 3) then  ! inverse DCT-3 is DCT-2
         lensav = 3*lenseq + 15
         allocate (wsave(lensav))
         call dcosqi(lenseq, wsave)
         call dcosqb(lenseq, result, wsave)
      end if
   end function idct_dp

end submodule fftpack_dct
