!*==RADF4.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RADF4(Ido,L1,Cc,Ch,Wa1,Wa2,Wa3)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RADF41932
!*** Start of declarations inserted by SPAG
      REAL Cc , Ch , ci2 , ci3 , ci4 , cr2 , cr3 , cr4 , FFTPACK_KIND , &
         & hsqt2 , rk , ti1 , ti2 , ti3 , ti4 , tr1 , tr2 , tr3 , tr4 , &
         & Wa1
      REAL Wa2 , Wa3
      INTEGER i , ic , Ido , idp2 , k , L1
!*** End of declarations inserted by SPAG
      DIMENSION Cc(Ido,L1,4) , Ch(Ido,4,L1) , Wa1(1) , Wa2(1) , Wa3(1)
      DATA hsqt2/0.70710678118654752440D0/
      DO k = 1 , L1
         tr1 = Cc(1,k,2) + Cc(1,k,4)
         tr2 = Cc(1,k,1) + Cc(1,k,3)
         Ch(1,1,k) = tr1 + tr2
         Ch(Ido,4,k) = tr2 - tr1
         Ch(Ido,2,k) = Cc(1,k,1) - Cc(1,k,3)
         Ch(1,3,k) = Cc(1,k,4) - Cc(1,k,2)
      ENDDO
      IF ( Ido<2 ) GOTO 99999
      IF ( Ido/=2 ) THEN
         idp2 = Ido + 2
         DO k = 1 , L1
            DO i = 3 , Ido , 2
               ic = idp2 - i
               cr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
               ci2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
               cr3 = Wa2(i-2)*Cc(i-1,k,3) + Wa2(i-1)*Cc(i,k,3)
               ci3 = Wa2(i-2)*Cc(i,k,3) - Wa2(i-1)*Cc(i-1,k,3)
               cr4 = Wa3(i-2)*Cc(i-1,k,4) + Wa3(i-1)*Cc(i,k,4)
               ci4 = Wa3(i-2)*Cc(i,k,4) - Wa3(i-1)*Cc(i-1,k,4)
               tr1 = cr2 + cr4
               tr4 = cr4 - cr2
               ti1 = ci2 + ci4
               ti4 = ci2 - ci4
               ti2 = Cc(i,k,1) + ci3
               ti3 = Cc(i,k,1) - ci3
               tr2 = Cc(i-1,k,1) + cr3
               tr3 = Cc(i-1,k,1) - cr3
               Ch(i-1,1,k) = tr1 + tr2
               Ch(ic-1,4,k) = tr2 - tr1
               Ch(i,1,k) = ti1 + ti2
               Ch(ic,4,k) = ti1 - ti2
               Ch(i-1,3,k) = ti4 + tr3
               Ch(ic-1,2,k) = tr3 - ti4
               Ch(i,3,k) = tr4 + ti3
               Ch(ic,2,k) = tr4 - ti3
            ENDDO
         ENDDO
         IF ( MOD(Ido,2)==1 ) RETURN
      ENDIF
      DO k = 1 , L1
         ti1 = -hsqt2*(Cc(Ido,k,2)+Cc(Ido,k,4))
         tr1 = hsqt2*(Cc(Ido,k,2)-Cc(Ido,k,4))
         Ch(Ido,1,k) = tr1 + Cc(Ido,k,1)
         Ch(Ido,3,k) = Cc(Ido,k,1) - tr1
         Ch(1,2,k) = ti1 - Cc(Ido,k,3)
         Ch(1,4,k) = ti1 + Cc(Ido,k,3)
      ENDDO
99999 END subroutine radf4