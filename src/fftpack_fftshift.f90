submodule(fftpack) fftpack_fftshift

contains

    !> Shifts zero-frequency component to center of spectrum for `complex` type.
    pure module function fftshift_crk(x) result(result)
        complex(kind=rk), intent(in) :: x(:)
        complex(kind=rk), allocatable :: result(:)

        result = cshift(x, shift=-floor(0.5_rk*size(x)))

    end function fftshift_crk

    !> Shifts zero-frequency component to center of spectrum for `real` type.
    pure module function fftshift_rrk(x) result(result)
        real(kind=rk), intent(in) :: x(:)
        real(kind=rk), allocatable :: result(:)

        result = cshift(x, shift=-floor(0.5_rk*size(x)))

    end function fftshift_rrk

end submodule fftpack_fftshift
