!*==DZFFTB.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DZFFTB(N,R,Azero,A,B,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DZFFTB567
!*** Start of declarations inserted by SPAG
      REAL A , Azero , B , FFTPACK_KIND , R , rk , Wsave
      INTEGER i , N , ns2
!*** End of declarations inserted by SPAG
      DIMENSION R(*) , A(*) , B(*) , Wsave(*)
      IF ( N<2 ) THEN
         R(1) = Azero
         RETURN
      ELSEIF ( N==2 ) THEN
         R(1) = Azero + A(1)
         R(2) = Azero - A(1)
         RETURN
      ELSE
         ns2 = (N-1)/2
         DO i = 1 , ns2
            R(2*i) = 0.5D0*A(i)
            R(2*i+1) = -0.5D0*B(i)
         ENDDO
         R(1) = Azero
         IF ( MOD(N,2)==0 ) R(N) = A(ns2+1)
         CALL DFFTB(N,R,Wsave(N+1))
      ENDIF
      END subroutine dzfftb