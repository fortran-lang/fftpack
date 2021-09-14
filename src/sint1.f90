!*==SINT1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE SINT1(N,War,Was,Xh,X,Ifac)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--SINT12466
!*** Start of declarations inserted by SPAG
      INTEGER i , Ifac , k , kc , modn , N , np1 , ns2
      REAL sqrt3 , t1 , t2 , War , Was , X , Xh , xhold
!*** End of declarations inserted by SPAG
      DIMENSION War(*) , Was(*) , X(*) , Xh(*) , Ifac(*)
      DATA sqrt3/1.73205080756887729352D0/
      DO i = 1 , N
         Xh(i) = War(i)
         War(i) = X(i)
      ENDDO
      IF ( N<2 ) THEN
         Xh(1) = Xh(1) + Xh(1)
      ELSEIF ( N==2 ) THEN
         xhold = sqrt3*(Xh(1)+Xh(2))
         Xh(2) = sqrt3*(Xh(1)-Xh(2))
         Xh(1) = xhold
      ELSE
         np1 = N + 1
         ns2 = N/2
         X(1) = 0.0D0
         DO k = 1 , ns2
            kc = np1 - k
            t1 = Xh(k) - Xh(kc)
            t2 = Was(k)*(Xh(k)+Xh(kc))
            X(k+1) = t1 + t2
            X(kc+1) = t2 - t1
         ENDDO
         modn = MOD(N,2)
         IF ( modn/=0 ) X(ns2+2) = 4.0D0*Xh(ns2+1)
         CALL RFFTF1(np1,X,Xh,War,Ifac)
         Xh(1) = 0.5D0*X(1)
         DO i = 3 , N , 2
            Xh(i-1) = -X(i)
            Xh(i) = Xh(i-2) + X(i-1)
         ENDDO
         IF ( modn==0 ) Xh(N) = -X(N+1)
      ENDIF
      DO i = 1 , N
         X(i) = War(i)
         War(i) = Xh(i)
      ENDDO
      END subroutine sint1