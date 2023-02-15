submodule(fftpack) fftpack_utils

contains

    !> Returns an integer array with the frequency values involved in the
    !> performed FFT, ordered in the standard way (zero first, then positive
    !> frequencies, then the negative ones).
    pure module function fftfreq(n) result(out)
        integer, intent(in) :: n
        integer, dimension(n) :: out
        integer :: i

        out(1) = 0
        if (n == 1) return

        if (mod(n, 2) == 0) then  !> n even, smallest n = 2
            do i = 2, n/2
                out(i) = i-1
            end do
            out(n/2+1) = -n/2
            do i = n/2+2, n  !> only enters if n/2+2 <= n
                out(i) = out(i-1) + 1
            end do
        else  !> n odd, smallest n = 3
            do i = 2, n/2+1
                out(i) = i-1
            end do
            out(n/2+2) = -out(n/2+1)
            do i = n/2+3, n  !> only enters if n/2+3 <= n
                out(i) = out(i-1) + 1
            end do
        end if
    end function fftfreq

    !> Returns an integer array with the frequency values involved in the
    !> performed real FFT, ordered in the standard way (zero first, then
    !> positive frequencies, then, if applicable, the negative one).
    pure module function rfftfreq(n) result(out)
        integer, intent(in) :: n
        integer, dimension(n) :: out
        integer :: i

        out(1) = 0
        if (n == 1) return

        if (mod(n,2) == 0) then  !> n even, smallest n = 2
            do i = 2, n-2, 2
                out(i) = out(i-1) + 1
                out(i+1) = out(i)
            end do
            out(n) = -n/2
        else  !> n odd, smallest n = 3
            do i = 2, n-1, 2
                out(i) = out(i-1) + 1
                out(i+1) = out(i)
            end do
        end if
    end function rfftfreq

end submodule fftpack_utils
