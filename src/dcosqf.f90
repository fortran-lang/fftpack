!*==DCOSQF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dcosqf(n,x,Wsave)
      use fftpack_kind
      implicit none
!*--DCOSQF311
!*** Start of declarations inserted by SPAG
      real fftpack_kind , rk , sqrt2 , tsqx , Wsave , x
      integer n
!*** End of declarations inserted by SPAG
      dimension x(*) , Wsave(*)
      data sqrt2/1.41421356237309504880d0/
      if ( n<2 ) then
      elseif ( n==2 ) then
         tsqx = sqrt2*x(2)
         x(2) = x(1) - tsqx
         x(1) = x(1) + tsqx
      else
         call cosqf1(n,x,Wsave,Wsave(n+1))
         goto 99999
      endif
      return
99999 end subroutine dcosqf