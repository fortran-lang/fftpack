!*==DSINQB.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dsinqb(n,x,Wsave)
      use fftpack_kind
      implicit none
!*--DSINQB469
!*** Start of declarations inserted by SPAG
      real fftpack_kind , rk , Wsave , x , xhold
      integer k , kc , n , ns2
!*** End of declarations inserted by SPAG
      dimension x(1) , Wsave(1)
      if ( n>1 ) then
         ns2 = n/2
         do k = 2 , n , 2
            x(k) = -x(k)
         enddo
         call dcosqb(n,x,Wsave)
         do k = 1 , ns2
            kc = n - k
            xhold = x(k)
            x(k) = x(kc+1)
            x(kc+1) = xhold
         enddo
         goto 99999
      endif
      x(1) = 4.0d0*x(1)
      return
99999 end subroutine dsinqb