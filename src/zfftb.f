      SUBROUTINE ZFFTB (N,C,WSAVE)
      USE fftpack_kind
      IMPLICIT REAL(RK) (A-H,O-Z)
      DIMENSION       C(1)       ,WSAVE(1)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
