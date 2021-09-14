!*==CFFTB1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE CFFTB1(N,C,Ch,Wa,Ifac)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--CFFTB15
!*** Start of declarations inserted by SPAG
      REAL C , Ch , FFTPACK_KIND , rk , Wa
      INTEGER i , idl1 , ido , idot , Ifac , ip , iw , ix2 , ix3 , ix4 ,&
            & k1 , l1 , l2 , N , n2 , na , nac , nf
!*** End of declarations inserted by SPAG
      DIMENSION Ch(*) , C(*) , Wa(*) , Ifac(*)
      nf = Ifac(2)
      na = 0
      l1 = 1
      iw = 1
      DO k1 = 1 , nf
         ip = Ifac(k1+2)
         l2 = ip*l1
         ido = N/l2
         idot = ido + ido
         idl1 = idot*l1
         IF ( ip==4 ) THEN
            ix2 = iw + idot
            ix3 = ix2 + idot
            IF ( na/=0 ) THEN
               CALL PASSB4(idot,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3))
            ELSE
               CALL PASSB4(idot,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3))
            ENDIF
            na = 1 - na
         ELSEIF ( ip==2 ) THEN
            IF ( na/=0 ) THEN
               CALL PASSB2(idot,l1,Ch,C,Wa(iw))
            ELSE
               CALL PASSB2(idot,l1,C,Ch,Wa(iw))
            ENDIF
            na = 1 - na
         ELSEIF ( ip==3 ) THEN
            ix2 = iw + idot
            IF ( na/=0 ) THEN
               CALL PASSB3(idot,l1,Ch,C,Wa(iw),Wa(ix2))
            ELSE
               CALL PASSB3(idot,l1,C,Ch,Wa(iw),Wa(ix2))
            ENDIF
            na = 1 - na
         ELSEIF ( ip/=5 ) THEN
            IF ( na/=0 ) THEN
               CALL PASSB(nac,idot,ip,l1,idl1,Ch,Ch,Ch,C,C,Wa(iw))
            ELSE
               CALL PASSB(nac,idot,ip,l1,idl1,C,C,C,Ch,Ch,Wa(iw))
            ENDIF
            IF ( nac/=0 ) na = 1 - na
         ELSE
            ix2 = iw + idot
            ix3 = ix2 + idot
            ix4 = ix3 + idot
            IF ( na/=0 ) THEN
               CALL PASSB5(idot,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
            ELSE
               CALL PASSB5(idot,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
            ENDIF
            na = 1 - na
         ENDIF
         l1 = l2
         iw = iw + (ip-1)*idot
      ENDDO
      IF ( na==0 ) RETURN
      n2 = N + N
      DO i = 1 , n2
         C(i) = Ch(i)
      ENDDO
      END subroutine cfftb1