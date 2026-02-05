      subroutine cosqf1(n, x, w, xh)
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
         do k = 2, ns2
            kc = np2 - k
            xh(k) = x(k) + x(kc)
            xh(kc) = x(k) - x(kc)
         end do
         modn = mod(n, 2)
         if (modn == 0) xh(ns2 + 1) = x(ns2 + 1) + x(ns2 + 1)
         do k = 2, ns2
            kc = np2 - k
            x(k) = w(k - 1)*xh(kc) + w(kc - 1)*xh(k)
            x(kc) = w(k - 1)*xh(k) - w(kc - 1)*xh(kc)
         end do
         if (modn == 0) x(ns2 + 1) = w(ns2)*xh(ns2 + 1)
         call dfftf(n, x, xh)
         do i = 3, n, 2
            xim1 = x(i - 1) - x(i)
            x(i) = x(i - 1) + x(i)
            x(i - 1) = xim1
         end do
      end subroutine cosqf1
