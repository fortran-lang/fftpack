      subroutine dcost(n, x, wsave)
         use fftpack_kind, only: rk
         implicit none
         integer, intent(in) :: n
         real(rk), intent(inout) :: wsave(*)
         real(rk), intent(inout) :: x(*)
         real(rk) :: c1, t1, t2, tx2, x1h, x1p3, &
                     xi, xim2
         integer :: i, k, kc, modn, nm1, np1, ns2
         nm1 = n - 1
         np1 = n + 1
         ns2 = n/2
         if (n < 2) return
         if (n == 2) then
            x1h = x(1) + x(2)
            x(2) = x(1) - x(2)
            x(1) = x1h
            return
         elseif (n > 3) then
            c1 = x(1) - x(n)
            x(1) = x(1) + x(n)
            do k = 2, ns2
               kc = np1 - k
               t1 = x(k) + x(kc)
               t2 = x(k) - x(kc)
               c1 = c1 + wsave(kc)*t2
               t2 = wsave(k)*t2
               x(k) = t1 - t2
               x(kc) = t1 + t2
            end do
            modn = mod(n, 2)
            if (modn /= 0) x(ns2 + 1) = x(ns2 + 1) + x(ns2 + 1)
            call dfftf(nm1, x, wsave(n + 1))
            xim2 = x(2)
            x(2) = c1
            do i = 4, n, 2
               xi = x(i)
               x(i) = x(i - 2) - x(i - 1)
               x(i - 1) = xim2
               xim2 = xi
            end do
            if (modn /= 0) x(n) = xim2
            return
         end if
         x1p3 = x(1) + x(3)
         tx2 = x(2) + x(2)
         x(2) = x(1) - x(3)
         x(1) = x1p3 + tx2
         x(3) = x1p3 - tx2
      end subroutine dcost
