!*==RFFTB1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RFFTB1(N,C,Ch,Wa,Ifac)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RFFTB12251
!*** Start of declarations inserted by SPAG
      REAL C , Ch , FFTPACK_KIND , rk , Wa
      INTEGER i , idl1 , ido , Ifac , ip , iw , ix2 , ix3 , ix4 , k1 ,  &
            & l1 , l2 , N , na , nf
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
         idl1 = ido*l1
         IF ( ip==4 ) THEN
            ix2 = iw + ido
            ix3 = ix2 + ido
            IF ( na/=0 ) THEN
               CALL RADB4(ido,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3))
            ELSE
               CALL RADB4(ido,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3))
            ENDIF
            na = 1 - na
         ELSEIF ( ip==2 ) THEN
            IF ( na/=0 ) THEN
               CALL RADB2(ido,l1,Ch,C,Wa(iw))
            ELSE
               CALL RADB2(ido,l1,C,Ch,Wa(iw))
            ENDIF
            na = 1 - na
         ELSEIF ( ip==3 ) THEN
            ix2 = iw + ido
            IF ( na/=0 ) THEN
               CALL RADB3(ido,l1,Ch,C,Wa(iw),Wa(ix2))
            ELSE
               CALL RADB3(ido,l1,C,Ch,Wa(iw),Wa(ix2))
            ENDIF
            na = 1 - na
         ELSEIF ( ip/=5 ) THEN
            IF ( na/=0 ) THEN
               CALL RADBG(ido,ip,l1,idl1,Ch,Ch,Ch,C,C,Wa(iw))
            ELSE
               CALL RADBG(ido,ip,l1,idl1,C,C,C,Ch,Ch,Wa(iw))
            ENDIF
            IF ( ido==1 ) na = 1 - na
         ELSE
            ix2 = iw + ido
            ix3 = ix2 + ido
            ix4 = ix3 + ido
            IF ( na/=0 ) THEN
               CALL RADB5(ido,l1,Ch,C,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
            ELSE
               CALL RADB5(ido,l1,C,Ch,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
            ENDIF
            na = 1 - na
         ENDIF
         l1 = l2
         iw = iw + (ip-1)*ido
      ENDDO
      IF ( na==0 ) RETURN
      DO i = 1 , N
         C(i) = Ch(i)
      ENDDO
      END subroutine rfftb1