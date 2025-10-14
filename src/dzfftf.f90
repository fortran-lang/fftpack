      subroutine dzfftf(n, r, azero, a, b, wsave)
!     version 3  june 1979
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(in) :: r(*)
         real(dp), intent(out) :: azero, a(*), b(*)
         real(dp), intent(inout) :: wsave(*)
         real(dp) :: cf, cfm
         integer :: i, ns2, ns2m
         if (n < 2) then
            azero = r(1)
            return
         elseif (n == 2) then
            azero = 0.5_dp*(r(1) + r(2))
            a(1) = 0.5_dp*(r(1) - r(2))
            return
         else
            do i = 1, n
               wsave(i) = r(i)
            end do
            call dfftf(n, wsave, wsave(n + 1))
            cf = 2.0_dp/real(n, dp)
            cfm = -cf
            azero = 0.5_dp*cf*wsave(1)
            ns2 = (n + 1)/2
            ns2m = ns2 - 1
            do i = 1, ns2m
               a(i) = cf*wsave(2*i)
               b(i) = cfm*wsave(2*i + 1)
            end do
            if (mod(n, 2) == 1) return
            a(ns2) = 0.5_dp*cf*wsave(n)
            b(ns2) = 0.0_dp
         end if
      end subroutine dzfftf
