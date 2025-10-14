      subroutine dzfftb(n, r, azero, a, b, wsave)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(out) :: r(*)
         real(dp), intent(inout) :: wsave(*)
         real(dp), intent(in) :: azero, a(*), b(*)
         integer :: i, ns2
         if (n < 2) then
            r(1) = azero
            return
         elseif (n == 2) then
            r(1) = azero + a(1)
            r(2) = azero - a(1)
            return
         else
            ns2 = (n - 1)/2
            do i = 1, ns2
               r(2*i) = 0.5_dp*a(i)
               r(2*i + 1) = -0.5_dp*b(i)
            end do
            r(1) = azero
            if (mod(n, 2) == 0) r(n) = a(ns2 + 1)
            call dfftb(n, r, wsave(n + 1))
         end if
      end subroutine dzfftb
