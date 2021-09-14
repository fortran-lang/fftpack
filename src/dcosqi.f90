!*==DCOSQI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DCOSQI(N,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DCOSQI333
!*** Start of declarations inserted by SPAG
      REAL dt , FFTPACK_KIND , fk , pih , rk , Wsave
      INTEGER k , N
!*** End of declarations inserted by SPAG
      DIMENSION Wsave(1)
      DATA pih/1.57079632679489661923D0/
      dt = pih/REAL(N,rk)
      fk = 0.0D0
      DO k = 1 , N
         fk = fk + 1.0D0
         Wsave(k) = COS(fk*dt)
      ENDDO
      CALL DFFTI(N,Wsave(N+1))
      END subroutine dcosqi