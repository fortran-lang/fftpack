!> fftpack kind
module fftpack_kind
    implicit none

    !> fftpack real kind
#if defined(fftpack_sp)
    integer, parameter :: rk = selected_real_kind(6)
#elif defined(fftpack_xdp)
    integer, parameter :: rk = selected_real_kind(18)
#elif defined(fftpack_qp)
    integer, parameter :: rk = selected_real_kind(33)
#else
    integer, parameter :: rk = selected_real_kind(15)
#endif

end module fftpack_kind
