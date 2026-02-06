      subroutine dsinti(n, wsave)
         use fftpack_kinds, only: dp
         implicit none
         integer, intent(in) :: n
         real(dp), intent(out) :: wsave(*)
         real(dp) :: dt
         integer :: k, np1, ns2
         real(dp), parameter :: pi = acos(-1.0_dp)
         if (n <= 1) return
         ns2 = n/2
         np1 = n + 1
         dt = pi/real(np1, dp)
         do k = 1, ns2
            wsave(k) = 2.0_dp*sin(k*dt)
         end do
         call dffti(np1, wsave(ns2 + 1))
      end subroutine dsinti
