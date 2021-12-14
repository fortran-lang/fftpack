submodule(fftpack) fftpack_dst

contains

    !> Discrete fourier sine transform of an odd sequence.
    pure module function dst_dp(x, n) result(result)
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
        lensav = int(2.5_dp*lenseq + 15)
        allocate (wsave(lensav))
        call dsinti(lenseq, wsave)

        !> Discrete fourier sine transformation
        call dsint(lenseq, result, wsave)

    end function dst_dp

end submodule fftpack_dst
