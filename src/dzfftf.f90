!*==DZFFTF.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dzfftf(n,r,Azero,a,b,Wsave)
!
!                       VERSION 3  JUNE 1979
!
      use fftpack_kind
      implicit none
!*--DZFFTF598
!*** Start of declarations inserted by SPAG
      real a , Azero , b , cf , cfm , fftpack_kind , r , rk , Wsave
      integer i , n , ns2 , ns2m
!*** End of declarations inserted by SPAG
      dimension r(*) , a(*) , b(*) , Wsave(*)
      if ( n<2 ) then
         Azero = r(1)
         return
      elseif ( n==2 ) then
         Azero = 0.5d0*(r(1)+r(2))
         a(1) = 0.5d0*(r(1)-r(2))
         return
      else
         do i = 1 , n
            Wsave(i) = r(i)
         enddo
         call dfftf(n,Wsave,Wsave(n+1))
         cf = 2.0d0/real(n,rk)
         cfm = -cf
         Azero = 0.5d0*cf*Wsave(1)
         ns2 = (n+1)/2
         ns2m = ns2 - 1
         do i = 1 , ns2m
            a(i) = cf*Wsave(2*i)
            b(i) = cfm*Wsave(2*i+1)
         enddo
         if ( mod(n,2)==1 ) return
         a(ns2) = 0.5d0*cf*Wsave(n)
         b(ns2) = 0.0d0
      endif
      end subroutine dzfftf