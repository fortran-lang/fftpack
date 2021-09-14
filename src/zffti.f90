!*==ZFFTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE ZFFTI(N,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--ZFFTI2543
!*** Start of declarations inserted by SPAG
      INTEGER iw1 , iw2 , N
      REAL Wsave
!*** End of declarations inserted by SPAG
      DIMENSION Wsave(1)
      IF ( N==1 ) RETURN
      iw1 = N + N + 1
      iw2 = iw1 + N + N
      CALL CFFTI1(N,Wsave(iw1),Wsave(iw2))
      END subroutine zffti