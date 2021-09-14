!*==DCOSTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DCOSTI(N,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DCOSTI405
!*** Start of declarations inserted by SPAG
      REAL dt , FFTPACK_KIND , fk , pi , rk , Wsave
      INTEGER k , kc , N , nm1 , np1 , ns2
!*** End of declarations inserted by SPAG
      DIMENSION Wsave(1)
      DATA pi/3.14159265358979323846D0/
      IF ( N<=3 ) RETURN
      nm1 = N - 1
      np1 = N + 1
      ns2 = N/2
      dt = pi/REAL(nm1,rk)
      fk = 0.0D0
      DO k = 2 , ns2
         kc = np1 - k
         fk = fk + 1.0D0
         Wsave(k) = 2.0D0*SIN(fk*dt)
         Wsave(kc) = 2.0D0*COS(fk*dt)
      ENDDO
      CALL DFFTI(nm1,Wsave(N+1))
      END subroutine dcosti