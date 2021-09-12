      SUBROUTINE DFFTB (N,R,WSAVE)
      USE fftpack_kind
      IMPLICIT REAL(RK) (A-H,O-Z)
      DIMENSION       R(1)       ,WSAVE(1)
      IF (N .EQ. 1) RETURN
      CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
