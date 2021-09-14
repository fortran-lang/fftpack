!*==RADF5.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RADF5(Ido,L1,Cc,Ch,Wa1,Wa2,Wa3,Wa4)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RADF51995
!*** Start of declarations inserted by SPAG
      REAL Cc , Ch , ci2 , ci3 , ci4 , ci5 , cr2 , cr3 , cr4 , cr5 ,    &
         & di2 , di3 , di4 , di5 , dr2 , dr3 , dr4 , dr5 ,              &
         & FFTPACK_KIND , rk
      REAL ti11 , ti12 , ti2 , ti3 , ti4 , ti5 , tr11 , tr12 , tr2 ,    &
         & tr3 , tr4 , tr5 , Wa1 , Wa2 , Wa3 , Wa4
      INTEGER i , ic , Ido , idp2 , k , L1
!*** End of declarations inserted by SPAG
      DIMENSION Cc(Ido,L1,5) , Ch(Ido,5,L1) , Wa1(1) , Wa2(1) , Wa3(1) ,&
              & Wa4(1)
      DATA tr11 , ti11 , tr12 , ti12/0.3090169943749474241D0 ,          &
         & 0.95105651629515357212D0 , -0.8090169943749474241D0 ,        &
         & 0.58778525229247312917D0/
      DO k = 1 , L1
         cr2 = Cc(1,k,5) + Cc(1,k,2)
         ci5 = Cc(1,k,5) - Cc(1,k,2)
         cr3 = Cc(1,k,4) + Cc(1,k,3)
         ci4 = Cc(1,k,4) - Cc(1,k,3)
         Ch(1,1,k) = Cc(1,k,1) + cr2 + cr3
         Ch(Ido,2,k) = Cc(1,k,1) + tr11*cr2 + tr12*cr3
         Ch(1,3,k) = ti11*ci5 + ti12*ci4
         Ch(Ido,4,k) = Cc(1,k,1) + tr12*cr2 + tr11*cr3
         Ch(1,5,k) = ti12*ci5 - ti11*ci4
      ENDDO
      IF ( Ido==1 ) RETURN
      idp2 = Ido + 2
      DO k = 1 , L1
         DO i = 3 , Ido , 2
            ic = idp2 - i
            dr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
            di2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
            dr3 = Wa2(i-2)*Cc(i-1,k,3) + Wa2(i-1)*Cc(i,k,3)
            di3 = Wa2(i-2)*Cc(i,k,3) - Wa2(i-1)*Cc(i-1,k,3)
            dr4 = Wa3(i-2)*Cc(i-1,k,4) + Wa3(i-1)*Cc(i,k,4)
            di4 = Wa3(i-2)*Cc(i,k,4) - Wa3(i-1)*Cc(i-1,k,4)
            dr5 = Wa4(i-2)*Cc(i-1,k,5) + Wa4(i-1)*Cc(i,k,5)
            di5 = Wa4(i-2)*Cc(i,k,5) - Wa4(i-1)*Cc(i-1,k,5)
            cr2 = dr2 + dr5
            ci5 = dr5 - dr2
            cr5 = di2 - di5
            ci2 = di2 + di5
            cr3 = dr3 + dr4
            ci4 = dr4 - dr3
            cr4 = di3 - di4
            ci3 = di3 + di4
            Ch(i-1,1,k) = Cc(i-1,k,1) + cr2 + cr3
            Ch(i,1,k) = Cc(i,k,1) + ci2 + ci3
            tr2 = Cc(i-1,k,1) + tr11*cr2 + tr12*cr3
            ti2 = Cc(i,k,1) + tr11*ci2 + tr12*ci3
            tr3 = Cc(i-1,k,1) + tr12*cr2 + tr11*cr3
            ti3 = Cc(i,k,1) + tr12*ci2 + tr11*ci3
            tr5 = ti11*cr5 + ti12*cr4
            ti5 = ti11*ci5 + ti12*ci4
            tr4 = ti12*cr5 - ti11*cr4
            ti4 = ti12*ci5 - ti11*ci4
            Ch(i-1,3,k) = tr2 + tr5
            Ch(ic-1,2,k) = tr2 - tr5
            Ch(i,3,k) = ti2 + ti5
            Ch(ic,2,k) = ti5 - ti2
            Ch(i-1,5,k) = tr3 + tr4
            Ch(ic-1,4,k) = tr3 - tr4
            Ch(i,5,k) = ti3 + ti4
            Ch(ic,4,k) = ti4 - ti3
         ENDDO
      ENDDO
      END subroutine radf5