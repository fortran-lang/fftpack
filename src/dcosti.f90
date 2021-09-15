!*==DCOSTI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dcosti(n,Wsave)
      use fftpack_kind
      implicit none
!*--DCOSTI405
!*** Start of declarations inserted by SPAG
      real dt , fftpack_kind , fk , pi , rk , Wsave
      integer k , kc , n , nm1 , np1 , ns2
!*** End of declarations inserted by SPAG
      dimension Wsave(1)
      data pi/3.14159265358979323846d0/
      if ( n<=3 ) return
      nm1 = n - 1
      np1 = n + 1
      ns2 = n/2
      dt = pi/real(nm1,rk)
      fk = 0.0d0
      do k = 2 , ns2
         kc = np1 - k
         fk = fk + 1.0d0
         Wsave(k) = 2.0d0*sin(fk*dt)
         Wsave(kc) = 2.0d0*cos(fk*dt)
      enddo
      call dffti(nm1,Wsave(n+1))
      end subroutine dcosti