!*==ZFFTF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine zfftf(n,c,Wsave)
      use fftpack_kind
      implicit none
!*--ZFFTF2528
!*** Start of declarations inserted by SPAG
      real c , Wsave
      integer iw1 , iw2 , n
!*** End of declarations inserted by SPAG
      dimension c(1) , Wsave(1)
      if ( n==1 ) return
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call cfftf1(n,c,Wsave,Wsave(iw1),Wsave(iw2))
      end subroutine zfftf