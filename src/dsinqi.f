      SUBROUTINE DSINQI (N,WSAVE)
      USE fftpack_kind
      IMPLICIT REAL(RK) (A-H,O-Z)
      DIMENSION       WSAVE(1)
      CALL DCOSQI (N,WSAVE)
      RETURN
      END
