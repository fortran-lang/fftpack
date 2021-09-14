!*==PASSF3.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE PASSF3(Ido,L1,Cc,Ch,Wa1,Wa2)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--PASSF31247
!*** Start of declarations inserted by SPAG
      REAL Cc , Ch , ci2 , ci3 , cr2 , cr3 , di2 , di3 , dr2 , dr3 ,    &
         & FFTPACK_KIND , rk , taui , taur , ti2 , tr2 , Wa1 , Wa2
      INTEGER i , Ido , k , L1
!*** End of declarations inserted by SPAG
      DIMENSION Cc(Ido,3,L1) , Ch(Ido,L1,3) , Wa1(1) , Wa2(1)
!     *** TAUI IS -SQRT(3)/2 ***
      DATA taur , taui/ - 0.5D0 , -0.86602540378443864676D0/
      IF ( Ido/=2 ) THEN
         DO k = 1 , L1
            DO i = 2 , Ido , 2
               tr2 = Cc(i-1,2,k) + Cc(i-1,3,k)
               cr2 = Cc(i-1,1,k) + taur*tr2
               Ch(i-1,k,1) = Cc(i-1,1,k) + tr2
               ti2 = Cc(i,2,k) + Cc(i,3,k)
               ci2 = Cc(i,1,k) + taur*ti2
               Ch(i,k,1) = Cc(i,1,k) + ti2
               cr3 = taui*(Cc(i-1,2,k)-Cc(i-1,3,k))
               ci3 = taui*(Cc(i,2,k)-Cc(i,3,k))
               dr2 = cr2 - ci3
               dr3 = cr2 + ci3
               di2 = ci2 + cr3
               di3 = ci2 - cr3
               Ch(i,k,2) = Wa1(i-1)*di2 - Wa1(i)*dr2
               Ch(i-1,k,2) = Wa1(i-1)*dr2 + Wa1(i)*di2
               Ch(i,k,3) = Wa2(i-1)*di3 - Wa2(i)*dr3
               Ch(i-1,k,3) = Wa2(i-1)*dr3 + Wa2(i)*di3
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
      DO k = 1 , L1
         tr2 = Cc(1,2,k) + Cc(1,3,k)
         cr2 = Cc(1,1,k) + taur*tr2
         Ch(1,k,1) = Cc(1,1,k) + tr2
         ti2 = Cc(2,2,k) + Cc(2,3,k)
         ci2 = Cc(2,1,k) + taur*ti2
         Ch(2,k,1) = Cc(2,1,k) + ti2
         cr3 = taui*(Cc(1,2,k)-Cc(1,3,k))
         ci3 = taui*(Cc(2,2,k)-Cc(2,3,k))
         Ch(1,k,2) = cr2 - ci3
         Ch(1,k,3) = cr2 + ci3
         Ch(2,k,2) = ci2 + cr3
         Ch(2,k,3) = ci2 - cr3
      ENDDO
      RETURN
99999 END subroutine passf3