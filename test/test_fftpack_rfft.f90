module test_fftpack_rfft

   use fftpack, only: dp, dffti, dfftf, dfftb, rfft, irfft
   use fftpack, only: dzffti, dzfftf, dzfftb
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_rfft

contains

   !> Collect all exported unit tests
   subroutine collect_rfft(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("classic-rfft-API", test_classic_rfft), &
                  new_unittest("modernized-rfft-API", test_modernized_rfft), &
                  new_unittest("modernized-irfft-API", test_modernized_irfft) &
                  ]

   end subroutine collect_rfft

   subroutine test_classic_rfft(error)
      type(error_type), allocatable, intent(out) :: error
      real(kind=dp) :: x(4) = [1, 2, 3, 4]
      real(kind=dp) :: w(31)
      real(kind=dp) :: azero, a(4/2), b(4/2)

      call dffti(4, w)
      call dfftf(4, x, w)
      call check(error, all(x == [real(kind=dp) :: 10, -2, 2, -2]), &
                 "`dfftf` failed.")
      if (allocated(error)) return
      call dfftb(4, x, w)
      call check(error, all(x/4.0_dp == [real(kind=dp) :: 1, 2, 3, 4]), &
                 "`dfftb` failed.")
      if (allocated(error)) return

      x = [1, 2, 3, 4]
      call dzffti(4, w)
      call dzfftf(4, x, azero, a, b, w)
      call check(error, azero == 2.5_dp, "dzfftf: azero == 2.5_dp failed.")
      if (allocated(error)) return
      call check(error, all(a == [-1.0_dp, -0.5_dp]), "dzfftf: all(a == [-1.0, -0.5]) failed.")
      if (allocated(error)) return
      call check(error, all(b == [-1.0_dp, 0.0_dp]), "dzfftf: all(b == [-1.0, 0.0]) failed.")
      if (allocated(error)) return

      call dzfftb(4, x, azero, a, b, w)
      call check(error, all(x == [real(kind=dp) :: 1, 2, 3, 4]), &
                 "dzfftb: all(x = [real(kind=dp) :: 1, 2, 3, 4]) failed.")

   end subroutine test_classic_rfft

   subroutine test_modernized_rfft(error)
      type(error_type), allocatable, intent(out) :: error
      real(kind=dp) :: eps = 1.0e-10_dp
      real(kind=dp) :: x(3) = [9, -9, 3]

      call check(error, sum(abs(rfft(x, 2) - [real(kind=dp) :: 0, 18])) < eps, &
                 "`rfft(x, 2)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(rfft(x, 3) - rfft(x))) < eps, &
                 "`rfft(x, 3)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(rfft(x, 4) - [real(kind=dp) :: 3, 6, 9, 21])) < eps, &
                 "`rfft(x, 4)` failed.")

   end subroutine test_modernized_rfft

   subroutine test_modernized_irfft(error)
      type(error_type), allocatable, intent(out) :: error
      real(kind=dp) :: eps = 1.0e-10_dp
      real(kind=dp) :: x(4) = [1, 2, 3, 4]

      call check(error, sum(abs(irfft(rfft(x))/4.0_dp - [real(kind=dp) :: 1, 2, 3, 4])) < eps, &
                 "`irfft(rfft(x))/4.0_dp` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(irfft(rfft(x), 2) - [real(kind=dp) :: 8, 12])) < eps, &
                 "`irfft(rfft(x), 2)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(irfft(rfft(x, 2), 4) - [real(kind=dp) :: 1, 3, 5, 3])) < eps, &
                 "`irfft(rfft(x, 2), 4)` failed.")

   end subroutine test_modernized_irfft

end module test_fftpack_rfft

program test_rfft
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_fftpack_rfft, only: collect_rfft
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("rfft", collect_rfft) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program
