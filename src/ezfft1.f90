!*==EZFFT1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE EZFFT1(N,Wa,Ifac)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--EZFFT1647
!*** Start of declarations inserted by SPAG
      REAL arg1 , argh , ch1 , ch1h , dch1 , dsh1 , FFTPACK_KIND , rk , &
         & sh1 , tpi , Wa
      INTEGER i , ib , ido , Ifac , ii , ip , ipm , is , j , k1 , l1 ,  &
            & l2 , N , nf , nfm1 , nl , nq , nr , ntry , ntryh
!*** End of declarations inserted by SPAG
      DIMENSION Wa(*) , Ifac(*) , ntryh(4)
      DATA ntryh(1) , ntryh(2) , ntryh(3) , ntryh(4)/4 , 2 , 3 , 5/ ,   &
         & tpi/6.28318530717958647692D0/
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
      argh = tpi/REAL(N,rk)
      is = 0
      nfm1 = nf - 1
      l1 = 1
      IF ( nfm1==0 ) RETURN
      DO k1 = 1 , nfm1
         ip = Ifac(k1+2)
         l2 = l1*ip
         ido = N/l2
         ipm = ip - 1
         arg1 = REAL(l1,rk)*argh
         ch1 = 1.0D0
         sh1 = 0.0D0
         dch1 = COS(arg1)
         dsh1 = SIN(arg1)
         DO j = 1 , ipm
            ch1h = dch1*ch1 - dsh1*sh1
            sh1 = dch1*sh1 + dsh1*ch1
            ch1 = ch1h
            i = is + 2
            Wa(i-1) = ch1
            Wa(i) = sh1
            IF ( ido>=5 ) THEN
               DO ii = 5 , ido , 2
                  i = i + 2
                  Wa(i-1) = ch1*Wa(i-3) - sh1*Wa(i-2)
                  Wa(i) = ch1*Wa(i-2) + sh1*Wa(i-3)
               ENDDO
            ENDIF
            is = is + ido
         ENDDO
         l1 = l2
      ENDDO
      END subroutine ezfft1