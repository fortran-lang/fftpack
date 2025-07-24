      subroutine dzfftf(n,r,Azero,a,b,Wsave)
!
!     VERSION 3  JUNE 1979
!
      use fftpack_kind
      implicit none
      real(rk) :: a , Azero , b , cf , cfm , r , Wsave
      integer :: i , n , ns2 , ns2m
      dimension r(*) , a(*) , b(*) , Wsave(*)
      if ( n<2 ) then
         Azero = r(1)
         return
      elseif ( n==2 ) then
         Azero = 0.5_rk*(r(1)+r(2))
         a(1) = 0.5_rk*(r(1)-r(2))
         return
      else
         do i = 1 , n
            Wsave(i) = r(i)
         enddo
         call dfftf(n,Wsave,Wsave(n+1))
         cf = 2.0_rk/real(n, rk)
         cfm = -cf
         Azero = 0.5_rk*cf*Wsave(1)
         ns2 = (n+1)/2
         ns2m = ns2 - 1
         do i = 1 , ns2m
            a(i) = cf*Wsave(2*i)
            b(i) = cfm*Wsave(2*i+1)
         enddo
         if ( mod(n,2)==1 ) return
         a(ns2) = 0.5_rk*cf*Wsave(n)
         b(ns2) = 0.0_rk
      endif
      end subroutine dzfftf