      subroutine dfftb(n, r, wsave)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(inout) :: r(*)
         real(dp), intent(inout) :: wsave(*)
         if (n == 1) return
         call rfftb1(n, r, wsave, wsave(n + 1), wsave(2*n + 1))
      end subroutine dfftb
