submodule(fftpack) fftpack_ifftshift

contains

    !> Shifts zero-frequency component to beginning of spectrum for `complex` type.
    pure module function ifftshift_crk(x) result(result)
        complex(kind=rk), intent(in) :: x(:)
        complex(kind=rk), allocatable :: result(:)

        result = cshift(x, shift=-ceiling(0.5_rk*size(x)))

    end function ifftshift_crk

    !> Shifts zero-frequency component to beginning of spectrum for `real` type.
    pure module function ifftshift_rrk(x) result(result)
        real(kind=rk), intent(in) :: x(:)
        real(kind=rk), allocatable :: result(:)

        result = cshift(x, shift=-ceiling(0.5_rk*size(x)))

    end function ifftshift_rrk

end submodule fftpack_ifftshift
