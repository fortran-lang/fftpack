!*==ZFFTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine zffti(n,Wsave)
      use fftpack_kind
      implicit none
!*--ZFFTI2543
!*** Start of declarations inserted by SPAG
      integer iw1 , iw2 , n
      real Wsave
!*** End of declarations inserted by SPAG
      dimension Wsave(1)
      if ( n==1 ) return
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call cffti1(n,Wsave(iw1),Wsave(iw2))
      end subroutine zffti