module test_fftpack_original
   use fftpack
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private

   public :: collect_original

   real(rk), parameter :: pi = 4.0_rk*atan(1.0_rk)
   real(rk), parameter :: atol = epsilon(1.0_rk)
   real(rk), parameter :: rtol = sqrt(atol)
   integer, parameter :: nd(1:7) = [120, 54, 49, 32, 4, 3, 2]

contains

   subroutine collect_original(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("dfft", test_dfft)]
      testsuite = [testsuite, new_unittest("zfft", test_zfft)]
      testsuite = [testsuite, new_unittest("sint", test_sint)]
   end subroutine collect_original

   subroutine test_dfft(error)
      type(error_type), allocatable, intent(out) :: error
      real(rk) :: x(200), y(200), xh(200), w(2000)
      integer :: i, j, k, n, np1, nm1, ns2, nz, modn
      real(rk) :: fn, tfn, dt, sum1, sum2, arg, arg1
      real(rk) :: mismatch, cf

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

         !> Check error.
         mismatch = maxval(abs(x(:n) - y(:n)))/fn
         call check(error, mismatch < rtol)
         if (allocated(error)) return

         !> Inverse Discrete Fourier Transform.
         x(:n) = xh(:n) ! Restore signal.
         do i = 1, n
            sum1 = 0.5_rk*x(1)
            arg = (i - 1)*dt
            if (ns2 >= 2) then
               do k = 2, ns2
                  arg1 = (k - 1)*arg
                  sum1 = sum1 + x(2*k - 2)*cos(arg1) - x(2*k - 1)*sin(arg1)
               end do
            end if
            if (modn == 0) sum1 = sum1 + 0.5_rk*(-1)**(i - 1)*x(n)
            y(i) = 2*sum1
         end do

         !> Perform (real) inverse Fourier Transform.
         call dfftb(n, x, w)

         !> Check error.
         mismatch = maxval(abs(x(:n) - y(:n)))
         call check(error, mismatch < rtol)
         if (allocated(error)) return

         !> Chain direct and inverse Fourier transforms.
         x(:n) = xh(:n); y(:n) = xh(:n) ! Restore signal.
         call dfftb(n, y, w)
         call dfftf(n, y, w)

         !> Check error.
         cf = 1.0_rk/fn
         mismatch = maxval(abs(cf*y(:n) - x(:n)))
         call check(error, mismatch < rtol)
         if (allocated(error)) return
      end do
   end subroutine test_dfft

   subroutine test_zfft(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j, k, n, nz
      complex(rk) :: cx(200), cy(200)
      real(rk) :: w(2000), dt, arg1, arg2, mismatch, cf

      do nz = 1, size(nd)
         !> Create signal.
         n = nd(nz)
         do i = 1, n
            cx(i) = cmplx(cos(sqrt(2.0_rk)*i), sin(sqrt(2.0_rk)*i**2), kind=rk)
         end do

         !> Discrete Fourier Transform.
         dt = 2*pi/n
         do i = 1, n
            arg1 = -(i - 1)*dt
            cy(i) = cmplx(0.0_rk, 0.0_rk, kind=rk)
            do k = 1, n
               arg2 = (k - 1)*arg1
               cy(i) = cy(i) + cmplx(cos(arg2), sin(arg2), kind=rk)*cx(k)
            end do
         end do

         !> Fast Fourier Transform.
         call zffti(n, w)
         call zfftf(n, cx, w)

         !> Check error.
         mismatch = maxval(abs(cx(:n) - cy(:n)))/n
         call check(error, mismatch < rtol)
         if (allocated(error)) return

         !> Inverse Discrete Fourier Transform.
         cx(:n) = cx(:n)/n ! Scale signal.
         do i = 1, n
            arg1 = (i - 1)*dt
            cy(i) = cmplx(0.0_rk, 0.0_rk, kind=rk)
            do k = 1, n
               arg2 = (k - 1)*arg1
               cy(i) = cy(i) + cmplx(cos(arg2), sin(arg2), kind=rk)*cx(k)
            end do
         end do

         !> Inverse Fast Fourier Transform.
         call zfftb(n, cx, w)

         !> Check error.
         mismatch = maxval(abs(cx(:n) - cy(:n)))
         call check(error, mismatch < rtol)
         if (allocated(error)) return

         !> Chain direct and inverse Fourier transforms.
         cx(:n) = cy(:n) ! Restore signal.
         call zfftf(n, cx, w)
         call zfftb(n, cx, w)

         !> Check error.
         cf = 1.0_rk/n
         mismatch = maxval(abs(cf*cx(:n) - cy(:n)))
         call check(error, mismatch < rtol)
         if (allocated(error)) return
      end do
   end subroutine test_zfft

   subroutine test_sint(error)
      type(error_type), allocatable, intent(out) :: error
      real(rk) :: x(200), y(200), xh(200), w(2000)
      integer :: i, j, k, n, np1, nm1, ns2, nz, modn
      real(rk) :: dt, sum1, sum2, arg, arg1
      real(rk) :: mismatch, cf

      do nz = 1, size(nd)
         !> Create multisine signal.
         n = nd(nz)
         modn = mod(n, 2)
         np1 = n + 1; nm1 = n - 1
         do j = 1, np1
            x(j) = sin(j*sqrt(2.0_rk))
            y(j) = x(j)
            xh(j) = x(j)
         end do

         !> Discrete sine transform.
         dt = pi/n
         do i = 1, nm1
            x(i) = xh(i)
         end do

         do i = 1, nm1
            y(i) = 0.0_rk
            arg1 = i*dt
            do k = 1, nm1
               y(i) = y(i) + x(k)*sin(k*arg1)
            end do
            y(i) = 2*y(i)
         end do

         !> Fast Sine Transform.
         call dsinti(nm1, w)
         call dsint(nm1, x, w)

         !> Check error.
         cf = 0.5_rk/n
         mismatch = maxval(abs(x(:nm1) - y(:nm1)))*cf
         call check(error, mismatch < rtol)
         if (allocated(error)) return

         !> Chain direct and inverse sine transform.
         x(:nm1) = xh(:nm1); y(:nm1) = x(:nm1) ! Restore signals.
         call dsint(nm1, x, w)
         call dsint(nm1, x, w)

         !> Check error.
         mismatch = maxval(abs(cf*x(:nm1) - y(:nm1)))
         call check(error, mismatch < rtol)
         if (allocated(error)) return
      end do
   end subroutine test_sint
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
