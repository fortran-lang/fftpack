      subroutine cosqb1(n, x, w, xh)
         use fftpack_kinds, only: dp
         implicit none
         integer, intent(in) :: n
         real(dp), intent(inout) :: x(*)
         real(dp), intent(in) :: w(*)
         real(dp), intent(out) :: xh(*)
         integer :: i, k, kc, modn, np2, ns2
         real(dp) :: xim1
         ns2 = (n + 1)/2
         np2 = n + 2
         do i = 3, n, 2
            xim1 = x(i - 1) + x(i)
            x(i) = x(i) - x(i - 1)
            x(i - 1) = xim1
         end do
         x(1) = x(1) + x(1)
         modn = mod(n, 2)
         if (modn == 0) x(n) = x(n) + x(n)
         call dfftb(n, x, xh)
         do k = 2, ns2
            kc = np2 - k
            xh(k) = w(k - 1)*x(kc) + w(kc - 1)*x(k)
            xh(kc) = w(k - 1)*x(k) - w(kc - 1)*x(kc)
         end do
         if (modn == 0) x(ns2 + 1) = w(ns2)*(x(ns2 + 1) + x(ns2 + 1))
         do k = 2, ns2
            kc = np2 - k
            x(k) = xh(k) + xh(kc)
            x(kc) = xh(k) - xh(kc)
         end do
         x(1) = x(1) + x(1)
      end subroutine cosqb1
