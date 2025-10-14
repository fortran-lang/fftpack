      subroutine dcosti(n, wsave)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(out) :: wsave(*)
         real(dp) :: dt, fk
         integer :: k, kc, nm1, np1, ns2
         real(dp), parameter :: pi = acos(-1.0_dp)
         if (n <= 3) return
         nm1 = n - 1
         np1 = n + 1
         ns2 = n/2
         dt = pi/real(nm1, kind=dp)
         fk = 0.0_dp
         do k = 2, ns2
            kc = np1 - k
            fk = fk + 1.0_dp
            wsave(k) = 2.0_dp*sin(fk*dt)
            wsave(kc) = 2.0_dp*cos(fk*dt)
         end do
         call dffti(nm1, wsave(n + 1))
      end subroutine dcosti
