      subroutine zfftb(n,c,Wsave)
      use fftpack_kind
      implicit none
      real(rk) :: c , Wsave
      integer :: iw1 , iw2 , n
      dimension c(1) , Wsave(*)
      if ( n==1 ) return
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call cfftb1(n,c,Wsave,Wsave(iw1),Wsave(iw2))
      end subroutine zfftb