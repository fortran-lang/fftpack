!*==COSQF1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE COSQF1(N,X,W,Xh)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--COSQF1256
!*** Start of declarations inserted by SPAG
      REAL FFTPACK_KIND , rk , W , X , Xh , xim1
      INTEGER i , k , kc , modn , N , np2 , ns2
!*** End of declarations inserted by SPAG
      DIMENSION X(1) , W(1) , Xh(1)
      ns2 = (N+1)/2
      np2 = N + 2
      DO k = 2 , ns2
         kc = np2 - k
         Xh(k) = X(k) + X(kc)
         Xh(kc) = X(k) - X(kc)
      ENDDO
      modn = MOD(N,2)
      IF ( modn==0 ) Xh(ns2+1) = X(ns2+1) + X(ns2+1)
      DO k = 2 , ns2
         kc = np2 - k
         X(k) = W(k-1)*Xh(kc) + W(kc-1)*Xh(k)
         X(kc) = W(k-1)*Xh(k) - W(kc-1)*Xh(kc)
      ENDDO
      IF ( modn==0 ) X(ns2+1) = W(ns2)*Xh(ns2+1)
      CALL DFFTF(N,X,Xh)
      DO i = 3 , N , 2
         xim1 = X(i-1) - X(i)
         X(i) = X(i-1) + X(i)
         X(i-1) = xim1
      ENDDO
      END subroutine cosqf1