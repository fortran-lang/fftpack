submodule(fftpack) fftpack_fft

contains

    !> Forward transform of a double complex periodic sequence.
    pure module function fft_cdp(x, n) result(result)
        complex(kind=dp), intent(in) :: x(:)
        integer, intent(in), optional :: n
        complex(kind=dp), allocatable :: result(:)

        integer :: lenseq, lensav, lenwrk
        real(kind=dp), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, ((0.0_dp, 0.0_dp), i=1, lenseq - size(x))]
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

    end function fft_cdp

end submodule fftpack_fft
