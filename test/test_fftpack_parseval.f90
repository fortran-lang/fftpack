module test_fftpack_parseval
   use fftpack, only: fft
   use fftpack_kind, only: rk
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_parseval
   integer, parameter :: n = 1024
   real(rk), parameter :: atol = n*epsilon(1.0_rk)
   real(rk), parameter :: rtol = sqrt(atol)

contains

   !> Collect all exported unit tests.
   subroutine collect_parseval(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Plancherel theorem", test_plancherel)]
      testsuite = [testsuite, new_unittest("Parseval theorem", test_parseval)]
   end subroutine collect_parseval

   subroutine test_plancherel(error)
      type(error_type), allocatable, intent(out) :: error
      !--------------------------------------
      !-----     COMPLEX(DP) TRANSFORM     -----
      !--------------------------------------
      block
         complex(rk) :: x(n)
         real(rk) :: xre(n), xim(n)
         real(rk) :: norms(2)
         call random_number(xre); call random_number(xim); x = cmplx(xre, xim, kind=rk)
         !> Norm of the signal in the time domain.
         norms(1) = sum(abs(x)**2)
         !> Normalize signal.
         x = x/sqrt(norms(1)); norms(1) = 1.0_rk
         !> Fourier transform.
         x = fft(x)
         !> Norm of the signal in frequency domain.
         norms(2) = sum(abs(x)**2)/n
         !> Check.
         call check(error, abs(norms(1) - norms(2)) < atol)
         if (allocated(error)) return
      end block
   end subroutine test_plancherel

   subroutine test_parseval(error)
      type(error_type), allocatable, intent(out) :: error
      !--------------------------------------
      !-----     COMPLEX(DP) TRANSFORM     -----
      !--------------------------------------
      block
         complex(rk) :: x(n), y(n), dotprods(2)
         real(rk) :: xre(n), xim(n)
         call random_number(xre); call random_number(xim); x = cmplx(xre, xim, kind=rk)
         call random_number(xre); call random_number(xim); y = cmplx(xre, xim, kind=rk)
         !> Normalize signals.
         x = x/sqrt(sum(abs(x)**2)); y = y/sqrt(sum(abs(y)**2))
         !> Dot product in time domain.
         dotprods(1) = dot_product(x, y)
         !> Fourier transform.
         x = fft(x); y = fft(y)
         !> Dot product in spectral domain.
         dotprods(2) = dot_product(x, y)/n
         !> Check.
         call check(error, abs(dotprods(1) - dotprods(2)) < atol)
         if (allocated(error)) return
      end block
   end subroutine test_parseval

end module test_fftpack_parseval

program test_parseval
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_fftpack_parseval, only: collect_parseval
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("fft", collect_parseval) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program test_parseval
