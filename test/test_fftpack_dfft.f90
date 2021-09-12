program tester

    call test_fftpack_dfft()
    print *, "All tests in `test_fftpack_dfft` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dfft()
        use fftpack, only: dffti, dfftf, dfftb
    use fftpack_kind

        real(kind=rk) :: x(4)
        real(kind=rk) :: w(31)

        x = [1, 2, 3, 4]
        
        call dffti(4, w)
        call dfftf(4, x, w)
        call check(all(x == [real(kind=rk) :: 10, -2, 2, -2]), &
                   msg="`dfftf` failed.")

        call dfftb(4, x, w)
        call check(all(x/4.0_rk == [real(kind=rk) :: 1, 2, 3, 4]), &
                   msg="`dfftb` failed.")

    end subroutine test_fftpack_dfft

end program tester
