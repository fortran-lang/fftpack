      SUBROUTINE DZFFTI (N,WSAVE)
      USE fftpack_kind
      IMPLICIT REAL(RK) (A-H,O-Z)
      DIMENSION       WSAVE(1)
      IF (N .EQ. 1) RETURN
      CALL EZFFT1 (N,WSAVE(2*N+1),WSAVE(3*N+1))
      RETURN
      END
