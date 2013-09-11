module types
implicit none
integer, parameter :: dp=kind(0.d0)
end module

!--------------------------------------------------------------------------

module constants
use types, only: dp
implicit none
real(dp), parameter :: pi    = 3.1415926535897932384626433832795_dp
end module

!--------------------------------------------------------------------------

module utils
use types, only: dp
implicit none

contains

subroutine init_random()
! Initializes the random number generator based on the system's time.
integer :: i, n, clock
integer, allocatable :: seed(:)
call random_seed(size=n)
allocate(seed(n))
call system_clock(count=clock)
seed = clock + 37 * [(i - 1, i = 1, n)]
call random_seed(put=seed)
end subroutine

subroutine assert(condition)
! If condition == .false., it aborts the program.
!
! Arguments
! ---------
!
logical, intent(in) :: condition
!
! Example
! -------
!
! call assert(a == 5)

if (.not. condition) call stop_error("Assert failed.")
end subroutine

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine

end module

!--------------------------------------------------------------------------

module fourier
use types, only: dp
use constants, only: pi
implicit none

contains

subroutine precalculate_coeffs(wa)
real(dp), intent(out) :: wa(:)
integer :: n, k, i, idx
n = size(wa) / 2
k = n / 2
wa = -1
wa(1) = 1
wa(2) = 0
idx = 2
do while (k > 0)
    do i = 1, k
        idx = idx + 2
        wa(idx-1) = cos(i*pi/k)
        wa(idx) = sin(i*pi/k)
        !print *, idx, i*pi/k, wa(idx-1), wa(idx)
    end do
    k = k/2
end do
end subroutine

subroutine passf2(IDO, L1, CC, CH, WA1)
integer, intent(in) :: IDO, L1
real(dp), intent(in) :: CC(IDO, 2, L1), WA1(:)
real(dp), intent(out) :: CH(IDO, L1, 2)
integer :: I, K
real(dp) :: TR2, TI2
if (IDO <= 2) then
    do K = 1, L1
        CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
        CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
        CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
        CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
    end do
else
    do K = 1, L1
        do I = 2, IDO, 2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
        end do
    end do
end if
end subroutine

end module

!--------------------------------------------------------------------------

program myf
use types, only: dp
use utils, only: init_random
use fourier, only: precalculate_coeffs
implicit none
real(dp) :: t1, t2, t3
real(dp), allocatable :: x_real(:), C(:)
complex(dp), allocatable :: x(:), xdft(:)
integer :: n
call init_random()
!n = 1024 * 1024
n = 32
allocate(x_real(n), x(n), xdft(n))
allocate(C(2*n))
call random_number(x_real)

x = x_real
call cpu_time(t1)
call precalculate_coeffs(C)
call cpu_time(t2)
call cpu_time(t3)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "zffti:", (t2-t1)*1000, "ms"
print *, "zfftf:", (t3-t2)*1000, "ms"
end program
