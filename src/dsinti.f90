!*==DSINTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dsinti(n,Wsave)
      use fftpack_kind
      implicit none
!*--DSINTI547
!*** Start of declarations inserted by SPAG
      real dt , fftpack_kind , pi , rk , Wsave
      integer k , n , np1 , ns2
!*** End of declarations inserted by SPAG
      dimension Wsave(1)
      data pi/3.14159265358979323846d0/
      if ( n<=1 ) return
      ns2 = n/2
      np1 = n + 1
      dt = pi/real(np1,rk)
      do k = 1 , ns2
         Wsave(k) = 2.0d0*sin(k*dt)
      enddo
      call dffti(np1,Wsave(ns2+1))
      end subroutine dsinti