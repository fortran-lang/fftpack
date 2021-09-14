!*==DCOST.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE DCOST(N,X,Wsave)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--DCOST352
!*** Start of declarations inserted by SPAG
      REAL c1 , FFTPACK_KIND , rk , t1 , t2 , tx2 , Wsave , X , x1h ,   &
         & x1p3 , xi , xim2
      INTEGER i , k , kc , modn , N , nm1 , np1 , ns2
!*** End of declarations inserted by SPAG
      DIMENSION X(*) , Wsave(*)
      nm1 = N - 1
      np1 = N + 1
      ns2 = N/2
      IF ( N<2 ) GOTO 99999
      IF ( N==2 ) THEN
         x1h = X(1) + X(2)
         X(2) = X(1) - X(2)
         X(1) = x1h
         RETURN
      ELSEIF ( N>3 ) THEN
         c1 = X(1) - X(N)
         X(1) = X(1) + X(N)
         DO k = 2 , ns2
            kc = np1 - k
            t1 = X(k) + X(kc)
            t2 = X(k) - X(kc)
            c1 = c1 + Wsave(kc)*t2
            t2 = Wsave(k)*t2
            X(k) = t1 - t2
            X(kc) = t1 + t2
         ENDDO
         modn = MOD(N,2)
         IF ( modn/=0 ) X(ns2+1) = X(ns2+1) + X(ns2+1)
         CALL DFFTF(nm1,X,Wsave(N+1))
         xim2 = X(2)
         X(2) = c1
         DO i = 4 , N , 2
            xi = X(i)
            X(i) = X(i-2) - X(i-1)
            X(i-1) = xim2
            xim2 = xi
         ENDDO
         IF ( modn/=0 ) X(N) = xim2
         GOTO 99999
      ENDIF
      x1p3 = X(1) + X(3)
      tx2 = X(2) + X(2)
      X(2) = X(1) - X(3)
      X(1) = x1p3 + tx2
      X(3) = x1p3 - tx2
      RETURN
99999 END subroutine dcost