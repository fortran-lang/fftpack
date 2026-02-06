      subroutine sint1(n, war, was, xh, x, ifac)
         use fftpack_kinds, only: dp
         implicit none
         integer, intent(in) :: n, ifac(*)
         real(dp), intent(in) :: was(*)
         real(dp), intent(inout) :: war(*), x(*)
         real(dp), intent(out) :: xh(*)
         integer :: i, k, kc, modn, np1, ns2
         real(dp) :: t1, t2, xhold
         real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
         do i = 1, n
            xh(i) = war(i)
            war(i) = x(i)
         end do
         if (n < 2) then
            xh(1) = xh(1) + xh(1)
         elseif (n == 2) then
            xhold = sqrt3*(xh(1) + xh(2))
            xh(2) = sqrt3*(xh(1) - xh(2))
            xh(1) = xhold
         else
            np1 = n + 1
            ns2 = n/2
            x(1) = 0.0_dp
            do k = 1, ns2
               kc = np1 - k
               t1 = xh(k) - xh(kc)
               t2 = was(k)*(xh(k) + xh(kc))
               x(k + 1) = t1 + t2
               x(kc + 1) = t2 - t1
            end do
            modn = mod(n, 2)
            if (modn /= 0) x(ns2 + 2) = 4.0_dp*xh(ns2 + 1)
            call rfftf1(np1, x, xh, war, ifac)
            xh(1) = 0.5_dp*x(1)
            do i = 3, n, 2
               xh(i - 1) = -x(i)
               xh(i) = xh(i - 2) + x(i - 1)
            end do
            if (modn == 0) xh(n) = -x(n + 1)
         end if
         do i = 1, n
            x(i) = war(i)
            war(i) = xh(i)
         end do
      end subroutine sint1
