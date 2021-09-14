!*==DZFFTF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DZFFTF(N,R,Azero,A,B,Wsave)
!
!                       VERSION 3  JUNE 1979
!
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DZFFTF598
!*** Start of declarations inserted by SPAG
      REAL A , Azero , B , cf , cfm , FFTPACK_KIND , R , rk , Wsave
      INTEGER i , N , ns2 , ns2m
!*** End of declarations inserted by SPAG
      DIMENSION R(*) , A(*) , B(*) , Wsave(*)
      IF ( N<2 ) THEN
         Azero = R(1)
         RETURN
      ELSEIF ( N==2 ) THEN
         Azero = 0.5D0*(R(1)+R(2))
         A(1) = 0.5D0*(R(1)-R(2))
         RETURN
      ELSE
         DO i = 1 , N
            Wsave(i) = R(i)
         ENDDO
         CALL DFFTF(N,Wsave,Wsave(N+1))
         cf = 2.0D0/REAL(N,rk)
         cfm = -cf
         Azero = 0.5D0*cf*Wsave(1)
         ns2 = (N+1)/2
         ns2m = ns2 - 1
         DO i = 1 , ns2m
            A(i) = cf*Wsave(2*i)
            B(i) = cfm*Wsave(2*i+1)
         ENDDO
         IF ( MOD(N,2)==1 ) RETURN
         A(ns2) = 0.5D0*cf*Wsave(N)
         B(ns2) = 0.0D0
      ENDIF
      END subroutine dzfftf