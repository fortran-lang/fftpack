!*==PASSF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE PASSF(Nac,Ido,Ip,L1,Idl1,Cc,C1,C2,Ch,Ch2,Wa)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--PASSF1087
!*** Start of declarations inserted by SPAG
      REAL C1 , C2 , Cc , Ch , Ch2 , FFTPACK_KIND , rk , Wa , wai , war
      INTEGER i , idij , idj , idl , Idl1 , idlj , Ido , idot , idp ,   &
            & ik , inc , Ip , ipp2 , ipph , j , jc , k , l , L1 , lc
      INTEGER Nac , nt
!*** End of declarations inserted by SPAG
      DIMENSION Ch(Ido,L1,Ip) , Cc(Ido,Ip,L1) , C1(Ido,L1,Ip) , Wa(1) , &
              & C2(Idl1,Ip) , Ch2(Idl1,Ip)
      idot = Ido/2
      nt = Ip*Idl1
      ipp2 = Ip + 2
      ipph = (Ip+1)/2
      idp = Ip*Ido
!
      IF ( Ido<L1 ) THEN
         DO j = 2 , ipph
            jc = ipp2 - j
            DO i = 1 , Ido
               DO k = 1 , L1
                  Ch(i,k,j) = Cc(i,j,k) + Cc(i,jc,k)
                  Ch(i,k,jc) = Cc(i,j,k) - Cc(i,jc,k)
               ENDDO
            ENDDO
         ENDDO
         DO i = 1 , Ido
            DO k = 1 , L1
               Ch(i,k,1) = Cc(i,1,k)
            ENDDO
         ENDDO
      ELSE
         DO j = 2 , ipph
            jc = ipp2 - j
            DO k = 1 , L1
               DO i = 1 , Ido
                  Ch(i,k,j) = Cc(i,j,k) + Cc(i,jc,k)
                  Ch(i,k,jc) = Cc(i,j,k) - Cc(i,jc,k)
               ENDDO
            ENDDO
         ENDDO
         DO k = 1 , L1
            DO i = 1 , Ido
               Ch(i,k,1) = Cc(i,1,k)
            ENDDO
         ENDDO
      ENDIF
      idl = 2 - Ido
      inc = 0
      DO l = 2 , ipph
         lc = ipp2 - l
         idl = idl + Ido
         DO ik = 1 , Idl1
            C2(ik,l) = Ch2(ik,1) + Wa(idl-1)*Ch2(ik,2)
            C2(ik,lc) = -Wa(idl)*Ch2(ik,Ip)
         ENDDO
         idlj = idl
         inc = inc + Ido
         DO j = 3 , ipph
            jc = ipp2 - j
            idlj = idlj + inc
            IF ( idlj>idp ) idlj = idlj - idp
            war = Wa(idlj-1)
            wai = Wa(idlj)
            DO ik = 1 , Idl1
               C2(ik,l) = C2(ik,l) + war*Ch2(ik,j)
               C2(ik,lc) = C2(ik,lc) - wai*Ch2(ik,jc)
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
         DO ik = 2 , Idl1 , 2
            Ch2(ik-1,j) = C2(ik-1,j) - C2(ik,jc)
            Ch2(ik-1,jc) = C2(ik-1,j) + C2(ik,jc)
            Ch2(ik,j) = C2(ik,j) + C2(ik-1,jc)
            Ch2(ik,jc) = C2(ik,j) - C2(ik-1,jc)
         ENDDO
      ENDDO
      Nac = 1
      IF ( Ido==2 ) RETURN
      Nac = 0
      DO ik = 1 , Idl1
         C2(ik,1) = Ch2(ik,1)
      ENDDO
      DO j = 2 , Ip
         DO k = 1 , L1
            C1(1,k,j) = Ch(1,k,j)
            C1(2,k,j) = Ch(2,k,j)
         ENDDO
      ENDDO
      IF ( idot>L1 ) THEN
         idj = 2 - Ido
         DO j = 2 , Ip
            idj = idj + Ido
            DO k = 1 , L1
               idij = idj
               DO i = 4 , Ido , 2
                  idij = idij + 2
                  C1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) + Wa(idij)       &
                              & *Ch(i,k,j)
                  C1(i,k,j) = Wa(idij-1)*Ch(i,k,j) - Wa(idij)           &
                            & *Ch(i-1,k,j)
               ENDDO
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
      idij = 0
      DO j = 2 , Ip
         idij = idij + 2
         DO i = 4 , Ido , 2
            idij = idij + 2
            DO k = 1 , L1
               C1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) + Wa(idij)*Ch(i,k,j)
               C1(i,k,j) = Wa(idij-1)*Ch(i,k,j) - Wa(idij)*Ch(i-1,k,j)
            ENDDO
         ENDDO
      ENDDO
      RETURN
99999 END subroutine passf