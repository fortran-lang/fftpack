      SUBROUTINE DSINT (N,X,WSAVE)
      USE fftpack_kind
      IMPLICIT REAL(RK) (A-H,O-Z)
      DIMENSION       X(1)       ,WSAVE(1)      
      NP1 = N+1
      IW1 = N/2+1
      IW2 = IW1+NP1
      IW3 = IW2+NP1
      CALL SINT1(N,X,WSAVE,WSAVE(IW1),WSAVE(IW2),WSAVE(IW3))
      RETURN
      END
