!*==CFFTI1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE CFFTI1(N,Wa,Ifac)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--CFFTI1149
!*** Start of declarations inserted by SPAG
      REAL arg , argh , argld , FFTPACK_KIND , fi , rk , tpi , Wa
      INTEGER i , i1 , ib , ido , idot , Ifac , ii , ip , ipm , j , k1 ,&
            & l1 , l2 , ld , N , nf , nl , nq , nr , ntry
      INTEGER ntryh
!*** End of declarations inserted by SPAG
      DIMENSION Wa(*) , Ifac(*) , ntryh(4)
      DATA ntryh(1) , ntryh(2) , ntryh(3) , ntryh(4)/3 , 4 , 2 , 5/
      nl = N
      nf = 0
      j = 0
 100  j = j + 1
      IF ( j<=4 ) THEN
         ntry = ntryh(j)
      ELSE
         ntry = ntry + 2
      ENDIF
 200  nq = nl/ntry
      nr = nl - ntry*nq
      IF ( nr/=0 ) GOTO 100
      nf = nf + 1
      Ifac(nf+2) = ntry
      nl = nq
      IF ( ntry==2 ) THEN
         IF ( nf/=1 ) THEN
            DO i = 2 , nf
               ib = nf - i + 2
               Ifac(ib+2) = Ifac(ib+1)
            ENDDO
            Ifac(3) = 2
         ENDIF
      ENDIF
      IF ( nl/=1 ) GOTO 200
      Ifac(1) = N
      Ifac(2) = nf
      tpi = 6.28318530717958647692D0
      argh = tpi/REAL(N,rk)
      i = 2
      l1 = 1
      DO k1 = 1 , nf
         ip = Ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = N/l2
         idot = ido + ido + 2
         ipm = ip - 1
         DO j = 1 , ipm
            i1 = i
            Wa(i-1) = 1.0D0
            Wa(i) = 0.0D0
            ld = ld + l1
            fi = 0.0D0
            argld = REAL(ld,rk)*argh
            DO ii = 4 , idot , 2
               i = i + 2
               fi = fi + 1.D0
               arg = fi*argld
               Wa(i-1) = COS(arg)
               Wa(i) = SIN(arg)
            ENDDO
            IF ( ip>5 ) THEN
               Wa(i1-1) = Wa(i-1)
               Wa(i1) = Wa(i)
            ENDIF
         ENDDO
         l1 = l2
      ENDDO
      END subroutine cffti1