program tester

    call test_fftpack_dzfft
    print *, "All tests in `test_fftpack_dzfft` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_dzfft
        use fftpack, only: dzffti, dzfftf, dzfftb
        use iso_fortran_env, only: dp => real64

        real(kind=dp) :: x(4) = [1, 2, 3, 4]
        real(kind=dp) :: w(3*4 + 15)
        real(kind=dp) :: azero, a(4/2), b(4/2)

        call dzffti(4, w)
        call dzfftf(4, x, azero, a, b, w)
        call check(azero == 2.5_dp, msg="azero == 2.5_dp failed.")
        call check(all(a == [real(kind=dp) :: -1.0, -0.5]), msg="all(a == [real(kind=dp) :: -1.0, -0.5]) failed.")
        call check(all(b == [real(kind=dp) :: -1.0, 0.0]), msg="all(b == [real(kind=dp) :: -1.0, 0.0]) failed.")
        
        x = 0
        call dzfftb(4, x, azero, a, b, w)
        call check(all(x == [real(kind=dp) :: 1, 2, 3, 4]), msg="all(x = [real(kind=dp) :: 1, 2, 3, 4]) failed.")

    end subroutine test_fftpack_dzfft

end program tester

