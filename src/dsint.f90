!*==DSINT.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dsint(n,x,Wsave)
      use fftpack_kind
      implicit none
!*--DSINT531
!*** Start of declarations inserted by SPAG
      real fftpack_kind , rk , Wsave , x
      integer iw1 , iw2 , iw3 , n , np1
!*** End of declarations inserted by SPAG
      dimension x(1) , Wsave(1)
      np1 = n + 1
      iw1 = n/2 + 1
      iw2 = iw1 + np1
      iw3 = iw2 + np1
      call sint1(n,x,Wsave,Wsave(iw1),Wsave(iw2),Wsave(iw3))
      end subroutine dsint