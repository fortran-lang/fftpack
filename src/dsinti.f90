!*==DSINTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DSINTI(N,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DSINTI547
!*** Start of declarations inserted by SPAG
      REAL dt , FFTPACK_KIND , pi , rk , Wsave
      INTEGER k , N , np1 , ns2
!*** End of declarations inserted by SPAG
      DIMENSION Wsave(1)
      DATA pi/3.14159265358979323846D0/
      IF ( N<=1 ) RETURN
      ns2 = N/2
      np1 = N + 1
      dt = pi/REAL(np1,rk)
      DO k = 1 , ns2
         Wsave(k) = 2.0D0*SIN(k*dt)
      ENDDO
      CALL DFFTI(np1,Wsave(ns2+1))
      END subroutine dsinti