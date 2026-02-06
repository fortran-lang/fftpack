      subroutine dsinqf(n, x, wsave)
         use fftpack_kinds, only: dp
         implicit none
         integer, intent(in) :: n
         real(dp), intent(inout) :: x(*)
         real(dp), intent(in) :: wsave(*)
         integer :: k, kc, ns2
         real(dp) :: xhold
         if (n == 1) return
         ns2 = n/2
         do k = 1, ns2
            kc = n - k
            xhold = x(k)
            x(k) = x(kc + 1)
            x(kc + 1) = xhold
         end do
         call dcosqf(n, x, wsave)
         do k = 2, n, 2
            x(k) = -x(k)
         end do
      end subroutine dsinqf
