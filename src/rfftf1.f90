!*==RFFTF1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RFFTF1(N,C,Ch,Wa,Ifac)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RFFTF12321
!*** Start of declarations inserted by SPAG
      REAL C , Ch , FFTPACK_KIND , rk , Wa
      INTEGER i , idl1 , ido , Ifac , ip , iw , ix2 , ix3 , ix4 , k1 ,  &
            & kh , l1 , l2 , N , na , nf
!*** End of declarations inserted by SPAG
      DIMENSION Ch(*) , C(*) , Wa(*) , Ifac(*)
      nf = Ifac(2)
      na = 1
      l2 = N
      iw = N
      DO k1 = 1 , nf
         kh = nf - k1
         ip = Ifac(kh+3)
         l1 = l2/ip
         ido = N/l2
         idl1 = ido*l1
         iw = iw - (ip-1)*ido
         na = 1 - na
         IF ( ip==4 ) THEN
            ix2 = iw + ido
            ix3 = ix2 + ido
            IF ( na/=0 ) THEN
               CALL RADF4(ido,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3))
            ELSE
               CALL RADF4(ido,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3))
            ENDIF
         ELSEIF ( ip/=2 ) THEN
            IF ( ip==3 ) THEN
               ix2 = iw + ido
               IF ( na/=0 ) THEN
                  CALL RADF3(ido,l1,Ch,C,Wa(iw),Wa(ix2))
               ELSE
                  CALL RADF3(ido,l1,C,Ch,Wa(iw),Wa(ix2))
               ENDIF
            ELSEIF ( ip/=5 ) THEN
               IF ( ido==1 ) na = 1 - na
               IF ( na/=0 ) THEN
                  CALL RADFG(ido,ip,l1,idl1,Ch,Ch,Ch,C,C,Wa(iw))
                  na = 0
               ELSE
                  CALL RADFG(ido,ip,l1,idl1,C,C,C,Ch,Ch,Wa(iw))
                  na = 1
               ENDIF
            ELSE
               ix2 = iw + ido
               ix3 = ix2 + ido
               ix4 = ix3 + ido
               IF ( na/=0 ) THEN
                  CALL RADF5(ido,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
               ELSE
                  CALL RADF5(ido,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
               ENDIF
            ENDIF
         ELSEIF ( na/=0 ) THEN
            CALL RADF2(ido,l1,Ch,C,Wa(iw))
         ELSE
            CALL RADF2(ido,l1,C,Ch,Wa(iw))
         ENDIF
         l2 = l1
      ENDDO
      IF ( na==1 ) RETURN
      DO i = 1 , N
         C(i) = Ch(i)
      ENDDO
      END subroutine rfftf1