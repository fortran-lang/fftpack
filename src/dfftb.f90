!*==DFFTB.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DFFTB(N,R,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DFFTB430
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , R , rk , Wsave
      INTEGER N
!*** End of declarations inserted by SPAG
      DIMENSION R(1) , Wsave(1)
      IF ( N==1 ) RETURN
      CALL RFFTB1(N,R,Wsave,Wsave(N+1),Wsave(2*N+1))
      END subroutine dfftb