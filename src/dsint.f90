      subroutine dsint(n, x, wsave)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(inout) :: x(*)
         real(dp), intent(in) :: wsave(*)
         integer :: iw1, iw2, iw3, np1
         np1 = n + 1
         iw1 = n/2 + 1
         iw2 = iw1 + np1
         iw3 = iw2 + np1
         call sint1(n, x, wsave, wsave(iw1), wsave(iw2), wsave(iw3))
      end subroutine dsint
