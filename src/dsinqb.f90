!*==DSINQB.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DSINQB(N,X,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DSINQB469
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , rk , Wsave , X , xhold
      INTEGER k , kc , N , ns2
!*** End of declarations inserted by SPAG
      DIMENSION X(1) , Wsave(1)
      IF ( N>1 ) THEN
         ns2 = N/2
         DO k = 2 , N , 2
            X(k) = -X(k)
         ENDDO
         CALL DCOSQB(N,X,Wsave)
         DO k = 1 , ns2
            kc = N - k
            xhold = X(k)
            X(k) = X(kc+1)
            X(kc+1) = xhold
         ENDDO
         GOTO 99999
      ENDIF
      X(1) = 4.0D0*X(1)
      RETURN
99999 END subroutine dsinqb