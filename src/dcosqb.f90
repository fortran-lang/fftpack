      subroutine dcosqb(n,x,Wsave)
      use fftpack_kind
      implicit none
      integer :: n
      real(rk) :: Wsave , x , x1
      dimension x(*) , Wsave(*)
      real(rk),parameter :: tsqrt2 = 2.0_rk * sqrt(2.0_rk)
      if ( n<2 ) then
         x(1) = 4.0_rk*x(1)
         return
      elseif ( n==2 ) then
         x1 = 4.0_rk*(x(1)+x(2))
         x(2) = tsqrt2*(x(1)-x(2))
         x(1) = x1
         return
      else
         call cosqb1(n,x,Wsave,Wsave(n+1))
      endif
      end subroutine dcosqb