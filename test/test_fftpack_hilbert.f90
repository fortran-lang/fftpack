program tester

    call test_fftpack_hilbert
    print *, "All tests in `test_fftpack_hilbert` passed."

contains

    subroutine check(condition, msg)
        logical, intent(in)          :: condition
        character(len=*), intent(in) :: msg
        if (.not. condition) error stop msg
    end subroutine check

    subroutine test_fftpack_hilbert()

        use fftpack, only: hilbert, dp
        real(kind=dp) :: eps = 1.0e-10_dp

        call check(sum(abs(hilbert([complex(kind=dp) :: 1, 2, 3, 4   ]) - &
                                   [complex(kind=dp) :: (1.0_dp, 1.0_dp),  (2.0_dp, -1.0_dp), &
                                                        (3.0_dp, -1.0_dp), (4.0_dp, 1.0_dp)])) < eps, &
                   msg="hilbert([complex(kind=dp) :: 1, 2, 3, 4]) failed.")
        call check(sum(abs(hilbert([complex(kind=dp) :: 1, 2, 3, 4, 5]) - &
                                   [complex(kind=dp) :: (1.0_dp, 1.7013016167040795_dp),   (2.0_dp, -1.3763819204711736_dp), &
                                                        (3.0_dp, -0.64983939246581257_dp), (4.0_dp, -1.3763819204711731_dp), &
                                                        (5.0_dp, 1.7013016167040800_dp)])) < eps, &
                   msg="hilbert([complex(kind=dp) :: 1, 2, 3, 4, 5]) failed.")

    end subroutine test_fftpack_hilbert

end program tester
