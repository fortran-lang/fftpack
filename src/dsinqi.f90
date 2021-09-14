      subroutine dsinqi(n,Wsave)
      use fftpack_kind
      implicit none
      integer :: n
      real(rk) :: Wsave
      dimension Wsave(*)
      call dcosqi(n,Wsave)
      end subroutine dsinqi