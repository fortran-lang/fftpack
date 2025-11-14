module test_fftpack_original
   use fftpack
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private

   public :: collect_original

   real(rk), parameter :: pi = 4.0_rk*atan(1.0_rk)
   real(rk), parameter :: atol = epsilon(1.0_rk)
   real(rk), parameter :: rtol = sqrt(atol)

contains

   subroutine collect_original(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("rfft", test_rfft)]
   end subroutine collect_original

   subroutine test_rfft(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nd(1:7) = [120, 54, 49, 32, 4, 3, 2]
      real(rk) :: x(200), y(200), xh(200), w(2000)
      integer :: i, j, k, n, np1, nm1, ns2, nz, modn
      real(rk) :: fn, tfn, dt, sum1, sum2, arg, arg1
      real(rk) :: mismatch

      do nz = 1, size(nd)
         !> Create multisine signal.
         n = nd(nz)
         modn = mod(n, 2)
         fn = real(n, kind=rk)
         tfn = 2*fn
         np1 = n + 1; nm1 = n - 1
         do j = 1, np1
            x(j) = sin(j*sqrt(2.0_rk))
            y(j) = x(j)
            xh(j) = x(j)
         end do

         !> Discrete Fourier Transform.
         dt = 2*pi/fn
         ns2 = (n + 1)/2
         if (ns2 >= 2) then
            do k = 2, ns2
               sum1 = 0.0_rk; sum2 = 0.0_rk
               arg = (k - 1)*dt
               do i = 1, n
                  arg1 = (i - 1)*arg
                  sum1 = sum1 + x(i)*cos(arg1)
                  sum2 = sum2 + x(i)*sin(arg1)
               end do
               y(2*k - 2) = sum1
               y(2*k - 1) = -sum2
            end do
         end if
         sum1 = 0.0_rk; sum2 = 0.0_rk
         do i = 1, nm1, 2
            sum1 = sum1 + x(i)
            sum2 = sum2 + x(i + 1)
         end do
         if (modn == 1) sum1 = sum1 + x(n)
         y(1) = sum1 + sum2
         if (modn == 0) y(n) = sum1 - sum2

         !> Perform (real) Fourier Transform.
         call dffti(n, w)
         call dfftf(n, x, w)

         !> Compute error.
         mismatch = 0.0_rk
         do i = 1, n
            mismatch = max(mismatch, abs(x(i) - y(i)))
            x(i) = xh(i)
         end do
         mismatch = mismatch/fn

         !> Check error.
         call check(error, mismatch < n*atol)
         if (allocated(error)) return
      end do
   end subroutine
end module test_fftpack_original

program test_original
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_fftpack_original, only: collect_original
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("Original test suite", collect_original) &
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
