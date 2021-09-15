!*==DCOSQB.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dcosqb(n,x,Wsave)
      use fftpack_kind
      implicit none
!*--DCOSQB288
!*** Start of declarations inserted by SPAG
      real fftpack_kind , rk , tsqrt2 , Wsave , x , x1
      integer n
!*** End of declarations inserted by SPAG
      dimension x(*) , Wsave(*)
      data tsqrt2/2.82842712474619009760d0/
      if ( n<2 ) then
         x(1) = 4.0d0*x(1)
         return
      elseif ( n==2 ) then
         x1 = 4.0d0*(x(1)+x(2))
         x(2) = tsqrt2*(x(1)-x(2))
         x(1) = x1
         return
      else
         call cosqb1(n,x,Wsave,Wsave(n+1))
      endif
      end subroutine dcosqb