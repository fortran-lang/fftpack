!*==DSINQI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DSINQI(N,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DSINQI519
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , rk , Wsave
      INTEGER N
!*** End of declarations inserted by SPAG
      DIMENSION Wsave(1)
      CALL DCOSQI(N,Wsave)
      END subroutine dsinqi