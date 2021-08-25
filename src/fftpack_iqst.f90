submodule(fftpack) fftpack_iqst

contains

    !> Backward sine-transform of quarter wave data.
    pure module function iqst_dp(x, n) result(result)
        real(kind=dp), intent(in) :: x(:)
        integer, intent(in), optional :: n
        real(kind=dp), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=dp), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, (0.0_dp, i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        !> Initialize FFT
        lensav = 3*lenseq + 15
        allocate (wsave(lensav))
        call dsinqi(lenseq, wsave)

        !> Backward transformation
        call dsinqb(lenseq, result, wsave)

    end function iqst_dp

end submodule fftpack_iqst
