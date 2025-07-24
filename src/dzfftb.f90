      subroutine dzfftb(n,r,Azero,a,b,Wsave)
      use fftpack_kind
      implicit none
      real(rk) :: a , Azero , b , r , Wsave
      integer :: i , n , ns2
      dimension r(*) , a(*) , b(*) , Wsave(*)
      if ( n<2 ) then
         r(1) = Azero
         return
      elseif ( n==2 ) then
         r(1) = Azero + a(1)
         r(2) = Azero - a(1)
         return
      else
         ns2 = (n-1)/2
         do i = 1 , ns2
            r(2*i) = 0.5_rk*a(i)
            r(2*i+1) = -0.5_rk*b(i)
         enddo
         r(1) = Azero
         if ( mod(n,2)==0 ) r(n) = a(ns2+1)
         call dfftb(n,r,Wsave(n+1))
      endif
      end subroutine dzfftb