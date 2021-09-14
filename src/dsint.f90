!*==DSINT.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DSINT(N,X,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DSINT531
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , rk , Wsave , X
      INTEGER iw1 , iw2 , iw3 , N , np1
!*** End of declarations inserted by SPAG
      DIMENSION X(1) , Wsave(1)
      np1 = N + 1
      iw1 = N/2 + 1
      iw2 = iw1 + np1
      iw3 = iw2 + np1
      CALL SINT1(N,X,Wsave,Wsave(iw1),Wsave(iw2),Wsave(iw3))
      END subroutine dsint