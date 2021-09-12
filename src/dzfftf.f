      SUBROUTINE DZFFTF (N,R,AZERO,A,B,WSAVE)
      USE fftpack_kind
C
C                       VERSION 3  JUNE 1979
C
      IMPLICIT REAL(RK) (A-H,O-Z)
      DIMENSION       R(*)       ,A(*)       ,B(*)       ,WSAVE(*)
      IF (N-2) 101,102,103
  101 AZERO = R(1)
      RETURN
  102 AZERO = 0.5D0*(R(1)+R(2))
      A(1) = 0.5D0*(R(1)-R(2))
      RETURN
  103 DO 104 I=1,N
         WSAVE(I) = R(I)
  104 CONTINUE
      CALL DFFTF (N,WSAVE,WSAVE(N+1))
      CF = 2.0D0/FLOAT(N)
      CFM = -CF
      AZERO = 0.5D0*CF*WSAVE(1)
      NS2 = (N+1)/2
      NS2M = NS2-1
      DO 105 I=1,NS2M
         A(I) = CF*WSAVE(2*I)
         B(I) = CFM*WSAVE(2*I+1)
  105 CONTINUE
      IF (MOD(N,2) .EQ. 1) RETURN
      A(NS2) = 0.5D0*CF*WSAVE(N)
      B(NS2) = 0.0D0
      RETURN
      END
