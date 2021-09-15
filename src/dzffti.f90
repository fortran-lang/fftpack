!*==DZFFTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dzffti(n,Wsave)
      use fftpack_kind
      implicit none
!*--DZFFTI634
!*** Start of declarations inserted by SPAG
      real fftpack_kind , rk , Wsave
      integer n
!*** End of declarations inserted by SPAG
      dimension Wsave(1)
      if ( n==1 ) return
      call ezfft1(n,Wsave(2*n+1),Wsave(3*n+1))
      end subroutine dzffti