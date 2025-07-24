      subroutine dsinti(n,Wsave)
      use fftpack_kind
      implicit none
      real(rk) :: dt , Wsave
      integer :: k , n , np1 , ns2
      dimension Wsave(*)
      real(rk),parameter :: pi = acos(-1.0_rk)
      if ( n<=1 ) return
      ns2 = n/2
      np1 = n + 1
      dt = pi/real(np1, rk)
      do k = 1 , ns2
         Wsave(k) = 2.0_rk*sin(k*dt)
      enddo
      call dffti(np1,Wsave(ns2+1))
      end subroutine dsinti