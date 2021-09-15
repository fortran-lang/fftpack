!> Reference: https://ww2.mathworks.cn/help/signal/ref/hilbert.html
submodule(fftpack) fftpack_hilbert

    use iso_fortran_env, only: int8

contains

    !> Computes the discrete-time analytic signal using Hilbert.
    pure module function hilbert_dp(x, n) result(result)
        complex(kind=dp), intent(in)  :: x(:)
        integer, intent(in), optional :: n
        complex(kind=dp), allocatable :: result(:)

        integer                         :: n_, lensav, i
        real(kind=dp), allocatable      :: wsave(:)
        integer(kind=int8), allocatable :: H(:)

        if (present(n)) then
            n_ = n
            if (n_ <= size(x)) then
                result = x(:n_)
            else if (n_ > size(x)) then
                result = [x, ((0.0_dp, 0.0_dp), i=1, n_ - size(x))]
            end if
        else
            n_ = size(x)
            result = x
        end if
        
        allocate(H(n_))

        !> Create a vector H whose elements H(i) have the values:
        !> 1 for i = 1
        !> 2 for i = 2, 3, ... , (n/2)
        !> 1(2) for i = n/2 + 1, mod(n, 2) ==(/=) 0
        !> 0 for i = (n/2)+2, ... , n
        H(1)           = 1_int8
        H(2:n_/2)      = 2_int8
        H(n_/2 + 1)    = merge(1_int8, 2_int8, mod(n_, 2) == 0)
        H(n_/2 + 2:n_) = 0_int8

        !> Initialize FFT
        lensav = 4*n_ + 15
        allocate(wsave(lensav))
        call zffti(n_, wsave)

        !> Forward transformation
        call zfftf(n_, result, wsave)
        
        result = result*H

        !> Backward transformation
        call zfftb(n_, result, wsave)

        result = result/n_

    end function hilbert_dp

end submodule fftpack_hilbert
