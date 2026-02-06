      subroutine dcosqb(n, x, wsave)
         use fftpack_kinds, only: dp
         implicit none
         integer, intent(in) :: n
         real(dp), intent(in) :: wsave(*)
         real(dp), intent(inout) :: x(*)
         real(dp) :: x1
         real(dp), parameter :: tsqrt2 = 2.0_dp*sqrt(2.0_dp)
         if (n < 2) then
            x(1) = 4.0_dp*x(1)
            return
         elseif (n == 2) then
            x1 = 4.0_dp*(x(1) + x(2))
            x(2) = tsqrt2*(x(1) - x(2))
            x(1) = x1
            return
         else
            call cosqb1(n, x, wsave, wsave(n + 1))
         end if
      end subroutine dcosqb
