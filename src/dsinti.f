      SUBROUTINE DSINTI (N,WSAVE)
      USE fftpack_kind
      IMPLICIT REAL(RK) (A-H,O-Z)
      DIMENSION       WSAVE(1)
      DATA PI /3.14159265358979323846D0/
      IF (N .LE. 1) RETURN
      NS2 = N/2
      NP1 = N+1
      DT = PI/REAL(NP1,RK)
      DO 101 K=1,NS2
         WSAVE(K) = 2.0D0*SIN(K*DT)
  101 CONTINUE
      CALL DFFTI (NP1,WSAVE(NS2+1))
      RETURN
      END
