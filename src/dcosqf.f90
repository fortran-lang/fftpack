!*==DCOSQF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DCOSQF(N,X,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DCOSQF311
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , rk , sqrt2 , tsqx , Wsave , X
      INTEGER N
!*** End of declarations inserted by SPAG
      DIMENSION X(*) , Wsave(*)
      DATA sqrt2/1.41421356237309504880D0/
      IF ( N<2 ) THEN
      ELSEIF ( N==2 ) THEN
         tsqx = sqrt2*X(2)
         X(2) = X(1) - tsqx
         X(1) = X(1) + tsqx
      ELSE
         CALL COSQF1(N,X,Wsave,Wsave(N+1))
         GOTO 99999
      ENDIF
      RETURN
99999 END subroutine dcosqf