!*==RFFTI1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RFFTI1(N,Wa,Ifac)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RFFTI12391
!*** Start of declarations inserted by SPAG
      REAL arg , argh , argld , FFTPACK_KIND , fi , rk , tpi , Wa
      INTEGER i , ib , ido , Ifac , ii , ip , ipm , is , j , k1 , l1 ,  &
            & l2 , ld , N , nf , nfm1 , nl , nq , nr , ntry
      INTEGER ntryh
!*** End of declarations inserted by SPAG
      DIMENSION Wa(*) , Ifac(*) , ntryh(4)
      DATA ntryh(1) , ntryh(2) , ntryh(3) , ntryh(4)/4 , 2 , 3 , 5/
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
      is = 0
      nfm1 = nf - 1
      l1 = 1
      IF ( nfm1==0 ) RETURN
      DO k1 = 1 , nfm1
         ip = Ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = N/l2
         ipm = ip - 1
         DO j = 1 , ipm
            ld = ld + l1
            i = is
            argld = REAL(ld,rk)*argh
            fi = 0.0D0
            DO ii = 3 , ido , 2
               i = i + 2
               fi = fi + 1.0D0
               arg = fi*argld
               Wa(i-1) = COS(arg)
               Wa(i) = SIN(arg)
            ENDDO
            is = is + ido
         ENDDO
         l1 = l2
      ENDDO
      END subroutine rffti1