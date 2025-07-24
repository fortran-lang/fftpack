module test_fftpack_dct

   use fftpack
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_dct

contains

   !> Collect all exported unit tests
   subroutine collect_dct(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("classic-dct-API", test_classic_dct), &
                  new_unittest("modernized-dct-API", test_modernized_dct), &
                  new_unittest("modernized-idct-API", test_modernized_idct) &
                  ]

   end subroutine collect_dct

   subroutine test_classic_dct(error)
      type(error_type), allocatable, intent(out) :: error
      real(kind=rk) :: w(3*4 + 15), w2(3*4 + 15)
      real(kind=rk) :: x(4) = [1, 2, 3, 4]
      real(kind=rk) :: x2(4)
      real(kind=rk) :: eps = 1.0e-10_rk

      x2 = x
      call dcosti(4, w)
      call dcost(4, x, w)
      call check(error, sum(abs(x - [real(kind=rk) :: 15, -4, 0, -1.0000000000000009_rk])) < eps, &
                 "`dcost` failed.")
      if (allocated(error)) return

      call dct_t1i(4, w2)
      call dct_t1(4, x2, w2)
      call check(error, maxval(abs(x2 - x)) < eps, "dct_t1 failed")
      if (allocated(error)) return

      call dcost(4, x, w)
      call check(error, sum(abs(x/(2*3) - [real(kind=rk) :: 1, 2, 3, 4])) < eps, &
                 "2nd `dcost` failed.")
      if (allocated(error)) return

      call dct_t1(4, x2, w2)
      call check(error, maxval(abs(x2 - x)) < eps, "2nd dct_t1 failed")
      if (allocated(error)) return

      x = [1, 2, 3, 4]
      x2 = x
      call dcosqi(4, w)
      call dcosqf(4, x, w)
      call check(error, sum(abs(x - [11.999626276085150_rk, -9.1029432177492193_rk, &
                                     2.6176618435106480_rk, -1.5143449018465791_rk])) < eps, &
                 "`dcosqf` failed.")
      if (allocated(error)) return

      call dct_t23i(4, w2)
      call dct_t3(4, x2, w2)
      call check(error, maxval(abs(x2 - x)) < eps, "dct_t3 failed")
      if (allocated(error)) return

      call dcosqb(4, x, w)
      call check(error, sum(abs(x/(4*4) - [real(kind=rk) :: 1, 2, 3, 4])) < eps, &
                 "`dcosqb` failed.")
      if (allocated(error)) return

      call dct_t2(4, x2, w2)
      call check(error, maxval(abs(x2 - x)) < eps, "dct_t2 failed")

   end subroutine test_classic_dct

   subroutine test_modernized_dct(error)
      type(error_type), allocatable, intent(out) :: error
      real(kind=rk) :: eps = 1.0e-10_rk
      real(kind=rk) :: x(3) = [9, -9, 3]

      ! DCT-1
      call check(error, sum(abs(dct(x, 2, 1) - [0.0_rk, 18.0_rk])) < eps, "`dct(x,2,1)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(dct(x, 3, 1) - dct(x, type=1))) < eps, "`dct(x,3,1)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(dct(x, 4, 1) - [real(kind=rk) :: -3, -3.0000000000000036_rk, 15, 33])) < eps, &
                 "`dct(x,4,1)` failed.")
      !DCT-2
      x = [9, -9, 3]
      call check(error, sum(abs(dct(x, 3, 2) - [12.0_rk, 20.784609690826525_rk, 60.0_rk])) < eps, &
                 "`dct(x,3,2)` failed.")
      call check(error, sum(abs(dct(x, 3) - [12.0_rk, 20.784609690826525_rk, 60.0_rk])) < eps, &
                 "`dct(x,3)` failed.")
      call check(error, sum(abs(dct(x, 4, 2) - [12.0_rk, 14.890858416882008_rk, 42.426406871192853_rk, &
                                                58.122821125684993_rk])) < eps, "`dct(x,4,2)` failed.")

      ! DCT-3
      x = [9, -9, 3]
      call check(error, sum(abs(dct(x, 2, 3) - [-3.7279220613578570_rk, 21.727922061357859_rk])) < eps, &
                 "`dct(x,2,3)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(dct(x, 3, 3) - dct(x, type=3))) < eps, &
                 "`dct(x,3,3)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(dct(x, 4, 3) - [-3.3871908980838743_rk, -2.1309424696909023_rk, &
                                                11.645661095452331_rk, 29.872472272322447_rk])) < eps, &
                 "`dct(x,4,3)` failed.")

   end subroutine test_modernized_dct

   subroutine test_modernized_idct(error)
      type(error_type), allocatable, intent(out) :: error
      real(kind=rk) :: eps = 1.0e-10_rk
      real(kind=rk) :: x(4) = [1, 2, 3, 4]

      ! inverse DCT-1
      call check(error, sum(abs(idct(dct(x, type=1), type=1)/(2*3) - x)) < eps, &
                 "`idct(dct(x,type=1),type=1)/(2*3)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(idct(dct(x, type=1), 2, 1)/(2*1) - [5.5_rk, 9.5_rk])) < eps, &
                 "`idct(dct(x,type=1), 2, 1)/(2*1)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(idct(dct(x, 2, 1), 4, 1)/(2*3) - [0.16666666666666666_rk, 0.33333333333333331_rk, &
                                                                  0.66666666666666663_rk, 0.83333333333333315_rk])) < eps, &
                 "`idct(dct(x,2,1), 4, 1)/(2*3)` failed.")

      ! inverse DCT-2
      x = [1, 2, 3, 4]
      call check(error, sum(abs(idct(dct(x, type=2))/(4*4) - x)) < eps, &
                 "`idct(dct(x, type=2))/(4*4)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(idct(dct(x, type=2), n=2) - [22.156460020898692_rk, 57.843539979101308_rk])) < eps, &
                 "`idct(dct(x, type=2),n=2)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(idct(dct(x, n=2, type=2), n=4) - [6.7737481404944937_rk, 9.8352155994152106_rk, &
                                                                  14.164784400584789_rk, 17.226251859505506_rk])) < eps, &
                 "`idct(dct(x, n=2, type=2), n=4)` failed.")

      ! inverse DCT-3
      x = [1, 2, 3, 4]
      call check(error, sum(abs(idct(dct(x, type=3), type=3)/(4*4) - x)) < eps, &
                 "`idct(dct(x, type=3), type=3)/(4*4)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(idct(dct(x, type=3), n=2, type=3)/(4*2) - &
                                [1.4483415291679655_rk, 7.4608849947753271_rk])) < eps, &
                 "`idct(dct(x, type=3), n=2, type=3)/(4*2)` failed.")
      if (allocated(error)) return
      call check(error, sum(abs(idct(dct(x, n=2, type=3), n=4, type=3)/(4*4) - &
                                [0.5_rk, 0.70932417358418376_rk, 1.0_rk, 0.78858050747473762_rk])) < eps, &
                 "`idct(dct(x, n=2, type=3), n=4, type=3)/(4*4)` failed.")

   end subroutine test_modernized_idct

end module test_fftpack_dct

program test_dct
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_fftpack_dct, only: collect_dct
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("linalg_norm", collect_dct) &
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
