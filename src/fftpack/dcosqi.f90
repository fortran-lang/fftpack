      subroutine dcosqi(n,Wsave)
      use fftpack_kind
      implicit none
      real(rk) :: dt , fk , Wsave
      integer :: k , n
      dimension Wsave(*)
      real(rk),parameter :: pih = acos(-1.0_rk) / 2.0_rk ! pi / 2
      dt = pih/real(n, rk)
      fk = 0.0_rk
      do k = 1 , n
         fk = fk + 1.0_rk
         Wsave(k) = cos(fk*dt)
      enddo
      call dffti(n,Wsave(n+1))
      end subroutine dcosqi