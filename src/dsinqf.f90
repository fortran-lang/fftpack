!*==DSINQF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DSINQF(N,X,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DSINQF496
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , rk , Wsave , X , xhold
      INTEGER k , kc , N , ns2
!*** End of declarations inserted by SPAG
      DIMENSION X(1) , Wsave(1)
      IF ( N==1 ) RETURN
      ns2 = N/2
      DO k = 1 , ns2
         kc = N - k
         xhold = X(k)
         X(k) = X(kc+1)
         X(kc+1) = xhold
      ENDDO
      CALL DCOSQF(N,X,Wsave)
      DO k = 2 , N , 2
         X(k) = -X(k)
      ENDDO
      END subroutine dsinqf