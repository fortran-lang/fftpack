      subroutine dcosqi(n, wsave)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(out) :: wsave(*)
         real(dp) :: dt, fk
         integer :: k
         real(dp), parameter :: pih = acos(-1.0_dp)/2.0_dp ! pi / 2
         dt = pih/real(n, kind=dp)
         fk = 0.0_dp
         do k = 1, n
            fk = fk + 1.0_dp
            wsave(k) = cos(fk*dt)
         end do
         call dffti(n, wsave(n + 1))
      end subroutine dcosqi
