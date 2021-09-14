      subroutine dcosqf(n,x,Wsave)
      use fftpack_kind
      implicit none
      integer :: n
      real(rk) :: tsqx , Wsave , x
      dimension x(*) , Wsave(*)
      real(rk),parameter :: sqrt2 = sqrt(2.0_rk)
      if ( n<2 ) then
         return
      elseif ( n==2 ) then
         tsqx = sqrt2*x(2)
         x(2) = x(1) - tsqx
         x(1) = x(1) + tsqx
      else
         call cosqf1(n,x,Wsave,Wsave(n+1))
      endif
      end subroutine dcosqf