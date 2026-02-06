      subroutine dsinqi(n, wsave)
         use fftpack_kinds, only: dp
         implicit none
         integer, intent(in) :: n
         real(dp), intent(out) :: wsave(*)
         call dcosqi(n, wsave)
      end subroutine dsinqi
