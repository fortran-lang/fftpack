      subroutine zffti(n, wsave)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(out) :: wsave(*)
         integer :: iw1, iw2
         if (n == 1) return
         iw1 = n + n + 1
         iw2 = iw1 + n + n
         call cffti1(n, wsave(iw1), wsave(iw2))
      end subroutine zffti
