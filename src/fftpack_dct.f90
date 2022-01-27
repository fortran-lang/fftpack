submodule(fftpack) fftpack_dct

contains

    !> Discrete fourier cosine transform of an even sequence.
    pure module function dct_rk(x, n) result(result)
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
        lensav = 3*lenseq + 15
        allocate (wsave(lensav))
        call dcosti(lenseq, wsave)

        !> Discrete fourier cosine transformation
        call dcost(lenseq, result, wsave)

    end function dct_rk

end submodule fftpack_dct
