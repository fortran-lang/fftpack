!*==DZFFTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DZFFTI(N,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DZFFTI634
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , rk , Wsave
      INTEGER N
!*** End of declarations inserted by SPAG
      DIMENSION Wsave(1)
      IF ( N==1 ) RETURN
      CALL EZFFT1(N,Wsave(2*N+1),Wsave(3*N+1))
      END subroutine dzffti