!*==RADBG.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RADBG(Ido,Ip,L1,Idl1,Cc,C1,C2,Ch,Ch2,Wa)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RADBG1676
!*** Start of declarations inserted by SPAG
      REAL ai1 , ai2 , ar1 , ar1h , ar2 , ar2h , arg , C1 , C2 , Cc ,   &
         & Ch , Ch2 , dc2 , dcp , ds2 , dsp , FFTPACK_KIND , rk , tpi , &
         & Wa
      INTEGER i , ic , idij , Idl1 , Ido , idp2 , ik , Ip , ipp2 ,      &
            & ipph , is , j , j2 , jc , k , l , L1 , lc , nbd
!*** End of declarations inserted by SPAG
      DIMENSION Ch(Ido,L1,Ip) , Cc(Ido,Ip,L1) , C1(Ido,L1,Ip) ,         &
              & C2(Idl1,Ip) , Ch2(Idl1,Ip) , Wa(1)
      DATA tpi/6.28318530717958647692D0/
      arg = tpi/REAL(Ip,rk)
      dcp = COS(arg)
      dsp = SIN(arg)
      idp2 = Ido + 2
      nbd = (Ido-1)/2
      ipp2 = Ip + 2
      ipph = (Ip+1)/2
      IF ( Ido<L1 ) THEN
         DO i = 1 , Ido
            DO k = 1 , L1
               Ch(i,k,1) = Cc(i,1,k)
            ENDDO
         ENDDO
      ELSE
         DO k = 1 , L1
            DO i = 1 , Ido
               Ch(i,k,1) = Cc(i,1,k)
            ENDDO
         ENDDO
      ENDIF
      DO j = 2 , ipph
         jc = ipp2 - j
         j2 = j + j
         DO k = 1 , L1
            Ch(1,k,j) = Cc(Ido,j2-2,k) + Cc(Ido,j2-2,k)
            Ch(1,k,jc) = Cc(1,j2-1,k) + Cc(1,j2-1,k)
         ENDDO
      ENDDO
      IF ( Ido/=1 ) THEN
         IF ( nbd<L1 ) THEN
            DO j = 2 , ipph
               jc = ipp2 - j
               DO i = 3 , Ido , 2
                  ic = idp2 - i
                  DO k = 1 , L1
                     Ch(i-1,k,j) = Cc(i-1,2*j-1,k) + Cc(ic-1,2*j-2,k)
                     Ch(i-1,k,jc) = Cc(i-1,2*j-1,k) - Cc(ic-1,2*j-2,k)
                     Ch(i,k,j) = Cc(i,2*j-1,k) - Cc(ic,2*j-2,k)
                     Ch(i,k,jc) = Cc(i,2*j-1,k) + Cc(ic,2*j-2,k)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO j = 2 , ipph
               jc = ipp2 - j
               DO k = 1 , L1
                  DO i = 3 , Ido , 2
                     ic = idp2 - i
                     Ch(i-1,k,j) = Cc(i-1,2*j-1,k) + Cc(ic-1,2*j-2,k)
                     Ch(i-1,k,jc) = Cc(i-1,2*j-1,k) - Cc(ic-1,2*j-2,k)
                     Ch(i,k,j) = Cc(i,2*j-1,k) - Cc(ic,2*j-2,k)
                     Ch(i,k,jc) = Cc(i,2*j-1,k) + Cc(ic,2*j-2,k)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      ar1 = 1.0D0
      ai1 = 0.0D0
      DO l = 2 , ipph
         lc = ipp2 - l
         ar1h = dcp*ar1 - dsp*ai1
         ai1 = dcp*ai1 + dsp*ar1
         ar1 = ar1h
         DO ik = 1 , Idl1
            C2(ik,l) = Ch2(ik,1) + ar1*Ch2(ik,2)
            C2(ik,lc) = ai1*Ch2(ik,Ip)
         ENDDO
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         DO j = 3 , ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h
            DO ik = 1 , Idl1
               C2(ik,l) = C2(ik,l) + ar2*Ch2(ik,j)
               C2(ik,lc) = C2(ik,lc) + ai2*Ch2(ik,jc)
            ENDDO
         ENDDO
      ENDDO
      DO j = 2 , ipph
         DO ik = 1 , Idl1
            Ch2(ik,1) = Ch2(ik,1) + Ch2(ik,j)
         ENDDO
      ENDDO
      DO j = 2 , ipph
         jc = ipp2 - j
         DO k = 1 , L1
            Ch(1,k,j) = C1(1,k,j) - C1(1,k,jc)
            Ch(1,k,jc) = C1(1,k,j) + C1(1,k,jc)
         ENDDO
      ENDDO
      IF ( Ido/=1 ) THEN
         IF ( nbd<L1 ) THEN
            DO j = 2 , ipph
               jc = ipp2 - j
               DO i = 3 , Ido , 2
                  DO k = 1 , L1
                     Ch(i-1,k,j) = C1(i-1,k,j) - C1(i,k,jc)
                     Ch(i-1,k,jc) = C1(i-1,k,j) + C1(i,k,jc)
                     Ch(i,k,j) = C1(i,k,j) + C1(i-1,k,jc)
                     Ch(i,k,jc) = C1(i,k,j) - C1(i-1,k,jc)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO j = 2 , ipph
               jc = ipp2 - j
               DO k = 1 , L1
                  DO i = 3 , Ido , 2
                     Ch(i-1,k,j) = C1(i-1,k,j) - C1(i,k,jc)
                     Ch(i-1,k,jc) = C1(i-1,k,j) + C1(i,k,jc)
                     Ch(i,k,j) = C1(i,k,j) + C1(i-1,k,jc)
                     Ch(i,k,jc) = C1(i,k,j) - C1(i-1,k,jc)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      IF ( Ido==1 ) RETURN
      DO ik = 1 , Idl1
         C2(ik,1) = Ch2(ik,1)
      ENDDO
      DO j = 2 , Ip
         DO k = 1 , L1
            C1(1,k,j) = Ch(1,k,j)
         ENDDO
      ENDDO
      IF ( nbd>L1 ) THEN
         is = -Ido
         DO j = 2 , Ip
            is = is + Ido
            DO k = 1 , L1
               idij = is
               DO i = 3 , Ido , 2
                  idij = idij + 2
                  C1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) - Wa(idij)       &
                              & *Ch(i,k,j)
                  C1(i,k,j) = Wa(idij-1)*Ch(i,k,j) + Wa(idij)           &
                            & *Ch(i-1,k,j)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         is = -Ido
         DO j = 2 , Ip
            is = is + Ido
            idij = is
            DO i = 3 , Ido , 2
               idij = idij + 2
               DO k = 1 , L1
                  C1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) - Wa(idij)       &
                              & *Ch(i,k,j)
                  C1(i,k,j) = Wa(idij-1)*Ch(i,k,j) + Wa(idij)           &
                            & *Ch(i-1,k,j)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      END subroutine radbg