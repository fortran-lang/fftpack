!*==DFFTF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dfftf(n,r,Wsave)
      use fftpack_kind
      implicit none
!*--DFFTF443
!*** Start of declarations inserted by SPAG
      real fftpack_kind , r , rk , Wsave
      integer n
!*** End of declarations inserted by SPAG
      dimension r(1) , Wsave(1)
      if ( n==1 ) return
      call rfftf1(n,r,Wsave,Wsave(n+1),Wsave(2*n+1))
      end subroutine dfftf