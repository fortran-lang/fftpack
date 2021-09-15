!*==DCOSQI.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine dcosqi(n,Wsave)
      use fftpack_kind
      implicit none
!*--DCOSQI333
!*** Start of declarations inserted by SPAG
      real dt , fftpack_kind , fk , pih , rk , Wsave
      integer k , n
!*** End of declarations inserted by SPAG
      dimension Wsave(1)
      data pih/1.57079632679489661923d0/
      dt = pih/real(n,rk)
      fk = 0.0d0
      do k = 1 , n
         fk = fk + 1.0d0
         Wsave(k) = cos(fk*dt)
      enddo
      call dffti(n,Wsave(n+1))
      end subroutine dcosqi