      subroutine dzffti(n,Wsave)
      use fftpack_kind
      implicit none
      integer :: n
      real(rk) :: Wsave
      dimension Wsave(*)
      if ( n==1 ) return
      call ezfft1(n,Wsave(2*n+1),Wsave(3*n+1))
      end subroutine dzffti