!*==RADF2.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RADF2(Ido,L1,Cc,Ch,Wa1)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RADF21854
!*** Start of declarations inserted by SPAG
      REAL Cc , Ch , FFTPACK_KIND , rk , ti2 , tr2 , Wa1
      INTEGER i , ic , Ido , idp2 , k , L1
!*** End of declarations inserted by SPAG
      DIMENSION Ch(Ido,2,L1) , Cc(Ido,L1,2) , Wa1(1)
      DO k = 1 , L1
         Ch(1,1,k) = Cc(1,k,1) + Cc(1,k,2)
         Ch(Ido,2,k) = Cc(1,k,1) - Cc(1,k,2)
      ENDDO
      IF ( Ido<2 ) GOTO 99999
      IF ( Ido/=2 ) THEN
         idp2 = Ido + 2
         DO k = 1 , L1
            DO i = 3 , Ido , 2
               ic = idp2 - i
               tr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
               ti2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
               Ch(i,1,k) = Cc(i,k,1) + ti2
               Ch(ic,2,k) = ti2 - Cc(i,k,1)
               Ch(i-1,1,k) = Cc(i-1,k,1) + tr2
               Ch(ic-1,2,k) = Cc(i-1,k,1) - tr2
            ENDDO
         ENDDO
         IF ( MOD(Ido,2)==1 ) RETURN
      ENDIF
      DO k = 1 , L1
         Ch(1,2,k) = -Cc(Ido,k,2)
         Ch(Ido,1,k) = Cc(Ido,k,1)
      ENDDO
99999 END subroutine radf2