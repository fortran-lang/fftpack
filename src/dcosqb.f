      SUBROUTINE DCOSQB (N,X,WSAVE)
      USE fftpack_kind
      IMPLICIT REAL(RK) (A-H,O-Z)
      DIMENSION       X(*)       ,WSAVE(*)
      DATA TSQRT2 /2.82842712474619009760D0/
      IF (N-2) 101,102,103
  101 X(1) = 4.0D0*X(1)
      RETURN
  102 X1 = 4.0D0*(X(1)+X(2))
      X(2) = TSQRT2*(X(1)-X(2))
      X(1) = X1
      RETURN
  103 CALL COSQB1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
