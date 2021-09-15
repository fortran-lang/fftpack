!*==DFFTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dffti(n,Wsave)
      use fftpack_kind
      implicit none
!*--DFFTI456
!*** Start of declarations inserted by SPAG
      real fftpack_kind , rk , Wsave
      integer n
!*** End of declarations inserted by SPAG
      dimension Wsave(1)
      if ( n==1 ) return
      call rffti1(n,Wsave(n+1),Wsave(2*n+1))
      end subroutine dffti