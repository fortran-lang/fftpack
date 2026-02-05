      subroutine dcosqf(n, x, wsave)
         use fftpack_kinds, only: dp
         implicit none
         integer, intent(in) :: n
         real(dp), intent(in) :: wsave(*)
         real(dp), intent(inout) :: x(*)
         real(dp) :: tsqx
         real(dp), parameter :: sqrt2 = sqrt(2.0_dp)
         if (n < 2) then
            return
         elseif (n == 2) then
            tsqx = sqrt2*x(2)
            x(2) = x(1) - tsqx
            x(1) = x(1) + tsqx
         else
            call cosqf1(n, x, wsave, wsave(n + 1))
         end if
      end subroutine dcosqf
