submodule(fftpack) fftpack_irfft
   implicit none
contains

   !> Backward transform of a real periodic sequence.
   pure module function irfft_dp(x, n) result(result)
      real(kind=dp), intent(in) :: x(:)
      integer, intent(in), optional :: n
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

      !> Initialize FFT
      lensav = 2*lenseq + 15
      allocate (wsave(lensav))
      call dffti(lenseq, wsave)

      !> Backward transformation
      call dfftb(lenseq, result, wsave)

   end function irfft_dp

end submodule fftpack_irfft
