submodule(fftpack) fftpack_fft

contains

    !> Forward transform of a complex periodic sequence.
    pure module function fft_rk(x, n) result(result)
        complex(kind=rk), intent(in) :: x(:)
        integer, intent(in), optional :: n
        complex(kind=rk), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=rk), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, ((0.0_rk, 0.0_rk), i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        !> Initialize FFT
        lensav = 4*lenseq + 15
        allocate (wsave(lensav))
        call zffti(lenseq, wsave)

        !> Forward transformation
        call zfftf(lenseq, result, wsave)

    end function fft_rk

end submodule fftpack_fft
