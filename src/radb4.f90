!*==RADB4.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE RADB4(Ido,L1,Cc,Ch,Wa1,Wa2,Wa3)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--RADB41532
!*** Start of declarations inserted by SPAG
      REAL Cc , Ch , ci2 , ci3 , ci4 , cr2 , cr3 , cr4 , FFTPACK_KIND , &
         & rk , sqrt2 , ti1 , ti2 , ti3 , ti4 , tr1 , tr2 , tr3 , tr4 , &
         & Wa1
      REAL Wa2 , Wa3
      INTEGER i , ic , Ido , idp2 , k , L1
!*** End of declarations inserted by SPAG
      DIMENSION Cc(Ido,4,L1) , Ch(Ido,L1,4) , Wa1(1) , Wa2(1) , Wa3(1)
      DATA sqrt2/1.41421356237309504880D0/
      DO k = 1 , L1
         tr1 = Cc(1,1,k) - Cc(Ido,4,k)
         tr2 = Cc(1,1,k) + Cc(Ido,4,k)
         tr3 = Cc(Ido,2,k) + Cc(Ido,2,k)
         tr4 = Cc(1,3,k) + Cc(1,3,k)
         Ch(1,k,1) = tr2 + tr3
         Ch(1,k,2) = tr1 - tr4
         Ch(1,k,3) = tr2 - tr3
         Ch(1,k,4) = tr1 + tr4
      ENDDO
      IF ( Ido<2 ) GOTO 99999
      IF ( Ido/=2 ) THEN
         idp2 = Ido + 2
         DO k = 1 , L1
            DO i = 3 , Ido , 2
               ic = idp2 - i
               ti1 = Cc(i,1,k) + Cc(ic,4,k)
               ti2 = Cc(i,1,k) - Cc(ic,4,k)
               ti3 = Cc(i,3,k) - Cc(ic,2,k)
               tr4 = Cc(i,3,k) + Cc(ic,2,k)
               tr1 = Cc(i-1,1,k) - Cc(ic-1,4,k)
               tr2 = Cc(i-1,1,k) + Cc(ic-1,4,k)
               ti4 = Cc(i-1,3,k) - Cc(ic-1,2,k)
               tr3 = Cc(i-1,3,k) + Cc(ic-1,2,k)
               Ch(i-1,k,1) = tr2 + tr3
               cr3 = tr2 - tr3
               Ch(i,k,1) = ti2 + ti3
               ci3 = ti2 - ti3
               cr2 = tr1 - tr4
               cr4 = tr1 + tr4
               ci2 = ti1 + ti4
               ci4 = ti1 - ti4
               Ch(i-1,k,2) = Wa1(i-2)*cr2 - Wa1(i-1)*ci2
               Ch(i,k,2) = Wa1(i-2)*ci2 + Wa1(i-1)*cr2
               Ch(i-1,k,3) = Wa2(i-2)*cr3 - Wa2(i-1)*ci3
               Ch(i,k,3) = Wa2(i-2)*ci3 + Wa2(i-1)*cr3
               Ch(i-1,k,4) = Wa3(i-2)*cr4 - Wa3(i-1)*ci4
               Ch(i,k,4) = Wa3(i-2)*ci4 + Wa3(i-1)*cr4
            ENDDO
         ENDDO
         IF ( MOD(Ido,2)==1 ) RETURN
      ENDIF
      DO k = 1 , L1
         ti1 = Cc(1,2,k) + Cc(1,4,k)
         ti2 = Cc(1,4,k) - Cc(1,2,k)
         tr1 = Cc(Ido,1,k) - Cc(Ido,3,k)
         tr2 = Cc(Ido,1,k) + Cc(Ido,3,k)
         Ch(Ido,k,1) = tr2 + tr2
         Ch(Ido,k,2) = sqrt2*(tr1-ti1)
         Ch(Ido,k,3) = ti2 + ti2
         Ch(Ido,k,4) = -sqrt2*(tr1+ti1)
      ENDDO
99999 END subroutine radb4