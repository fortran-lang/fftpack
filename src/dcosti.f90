      subroutine dcosti(n,Wsave)
      use fftpack_kind
      implicit none
      real(rk) :: dt , fk , Wsave
      integer :: k , kc , n , nm1 , np1 , ns2
      dimension Wsave(*)
      real(rk),parameter :: pi = acos(-1.0_rk)
      if ( n<=3 ) return
      nm1 = n - 1
      np1 = n + 1
      ns2 = n/2
      dt = pi/real(nm1, rk)
      fk = 0.0_rk
      do k = 2 , ns2
         kc = np1 - k
         fk = fk + 1.0_rk
         Wsave(k) = 2.0_rk*sin(fk*dt)
         Wsave(kc) = 2.0_rk*cos(fk*dt)
      enddo
      call dffti(nm1,Wsave(n+1))
      end subroutine dcosti