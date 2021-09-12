submodule(fftpack) fftpack_irfft

contains

    !> Backward transform of a real periodic sequence.
    pure module function irfft_rk(x, n) result(result)
        real(kind=rk), intent(in) :: x(:)
        integer, intent(in), optional :: n
        real(kind=rk), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=rk), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, (0.0_rk, i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        !> Initialize FFT
        lensav = 2*lenseq + 15
        allocate (wsave(lensav))
        call dffti(lenseq, wsave)

        !> Backward transformation
        call dfftb(lenseq, result, wsave)

    end function irfft_rk

end submodule fftpack_irfft
