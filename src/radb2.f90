!*==RADB2.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RADB2(Ido,L1,Cc,Ch,Wa1)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RADB21452
!*** Start of declarations inserted by SPAG
      REAL Cc , Ch , FFTPACK_KIND , rk , ti2 , tr2 , Wa1
      INTEGER i , ic , Ido , idp2 , k , L1
!*** End of declarations inserted by SPAG
      DIMENSION Cc(Ido,2,L1) , Ch(Ido,L1,2) , Wa1(1)
      DO k = 1 , L1
         Ch(1,k,1) = Cc(1,1,k) + Cc(Ido,2,k)
         Ch(1,k,2) = Cc(1,1,k) - Cc(Ido,2,k)
      ENDDO
      IF ( Ido<2 ) GOTO 99999
      IF ( Ido/=2 ) THEN
         idp2 = Ido + 2
         DO k = 1 , L1
            DO i = 3 , Ido , 2
               ic = idp2 - i
               Ch(i-1,k,1) = Cc(i-1,1,k) + Cc(ic-1,2,k)
               tr2 = Cc(i-1,1,k) - Cc(ic-1,2,k)
               Ch(i,k,1) = Cc(i,1,k) - Cc(ic,2,k)
               ti2 = Cc(i,1,k) + Cc(ic,2,k)
               Ch(i-1,k,2) = Wa1(i-2)*tr2 - Wa1(i-1)*ti2
               Ch(i,k,2) = Wa1(i-2)*ti2 + Wa1(i-1)*tr2
            ENDDO
         ENDDO
         IF ( MOD(Ido,2)==1 ) RETURN
      ENDIF
      DO k = 1 , L1
         Ch(Ido,k,1) = Cc(Ido,1,k) + Cc(Ido,1,k)
         Ch(Ido,k,2) = -(Cc(1,2,k)+Cc(1,2,k))
      ENDDO
99999 END subroutine radb2