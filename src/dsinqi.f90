!*==DSINQI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dsinqi(n,Wsave)
      use fftpack_kind
      implicit none
!*--DSINQI519
!*** Start of declarations inserted by SPAG
      real fftpack_kind , rk , Wsave
      integer n
!*** End of declarations inserted by SPAG
      dimension Wsave(1)
      call dcosqi(n,Wsave)
      end subroutine dsinqi