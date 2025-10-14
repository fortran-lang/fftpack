      subroutine dsinqi(n, wsave)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(out) :: wsave(*)
         call dcosqi(n, wsave)
      end subroutine dsinqi
