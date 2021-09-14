!*==ZFFTF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE ZFFTF(N,C,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--ZFFTF2528
!*** Start of declarations inserted by SPAG
      REAL C , Wsave
      INTEGER iw1 , iw2 , N
!*** End of declarations inserted by SPAG
      DIMENSION C(1) , Wsave(1)
      IF ( N==1 ) RETURN
      iw1 = N + N + 1
      iw2 = iw1 + N + N
      CALL CFFTF1(N,C,Wsave,Wsave(iw1),Wsave(iw2))
      END subroutine zfftf