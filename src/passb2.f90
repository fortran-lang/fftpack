!*==PASSB2.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      SUBROUTINE PASSB2(Ido,L1,Cc,Ch,Wa1)
      USE FFTPACK_KIND
      IMPLICIT NONE
!*--PASSB2851
!*** Start of declarations inserted by SPAG
      REAL Cc , Ch , FFTPACK_KIND , rk , ti2 , tr2 , Wa1
      INTEGER i , Ido , k , L1
!*** End of declarations inserted by SPAG
      DIMENSION Cc(Ido,2,L1) , Ch(Ido,L1,2) , Wa1(1)
      IF ( Ido>2 ) THEN
         DO k = 1 , L1
            DO i = 2 , Ido , 2
               Ch(i-1,k,1) = Cc(i-1,1,k) + Cc(i-1,2,k)
               tr2 = Cc(i-1,1,k) - Cc(i-1,2,k)
               Ch(i,k,1) = Cc(i,1,k) + Cc(i,2,k)
               ti2 = Cc(i,1,k) - Cc(i,2,k)
               Ch(i,k,2) = Wa1(i-1)*ti2 + Wa1(i)*tr2
               Ch(i-1,k,2) = Wa1(i-1)*tr2 - Wa1(i)*ti2
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
      DO k = 1 , L1
         Ch(1,k,1) = Cc(1,1,k) + Cc(1,2,k)
         Ch(1,k,2) = Cc(1,1,k) - Cc(1,2,k)
         Ch(2,k,1) = Cc(2,1,k) + Cc(2,2,k)
         Ch(2,k,2) = Cc(2,1,k) - Cc(2,2,k)
      ENDDO
      RETURN
99999 END subroutine passb2