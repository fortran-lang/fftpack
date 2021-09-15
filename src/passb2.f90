!*==PASSB2.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine passb2(Ido,l1,Cc,Ch,Wa1)
      use fftpack_kind
      implicit none
!*--PASSB2851
!*** Start of declarations inserted by SPAG
      real Cc , Ch , fftpack_kind , rk , ti2 , tr2 , Wa1
      integer i , Ido , k , l1
!*** End of declarations inserted by SPAG
      dimension Cc(Ido,2,l1) , Ch(Ido,l1,2) , Wa1(1)
      if ( Ido>2 ) then
         do k = 1 , l1
            do i = 2 , Ido , 2
               Ch(i-1,k,1) = Cc(i-1,1,k) + Cc(i-1,2,k)
               tr2 = Cc(i-1,1,k) - Cc(i-1,2,k)
               Ch(i,k,1) = Cc(i,1,k) + Cc(i,2,k)
               ti2 = Cc(i,1,k) - Cc(i,2,k)
               Ch(i,k,2) = Wa1(i-1)*ti2 + Wa1(i)*tr2
               Ch(i-1,k,2) = Wa1(i-1)*tr2 - Wa1(i)*ti2
            enddo
         enddo
         goto 99999
      endif
      do k = 1 , l1
         Ch(1,k,1) = Cc(1,1,k) + Cc(1,2,k)
         Ch(1,k,2) = Cc(1,1,k) - Cc(1,2,k)
         Ch(2,k,1) = Cc(2,1,k) + Cc(2,2,k)
         Ch(2,k,2) = Cc(2,1,k) - Cc(2,2,k)
      enddo
      return
99999 end subroutine passb2