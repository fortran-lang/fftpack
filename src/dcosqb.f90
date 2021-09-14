!*==DCOSQB.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DCOSQB(N,X,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DCOSQB288
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , rk , tsqrt2 , Wsave , X , x1
      INTEGER N
!*** End of declarations inserted by SPAG
      DIMENSION X(*) , Wsave(*)
      DATA tsqrt2/2.82842712474619009760D0/
      IF ( N<2 ) THEN
         X(1) = 4.0D0*X(1)
         RETURN
      ELSEIF ( N==2 ) THEN
         x1 = 4.0D0*(X(1)+X(2))
         X(2) = tsqrt2*(X(1)-X(2))
         X(1) = x1
         RETURN
      ELSE
         CALL COSQB1(N,X,Wsave,Wsave(N+1))
      ENDIF
      END subroutine dcosqb