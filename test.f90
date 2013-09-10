program test
integer, parameter :: dp=kind(0.d0)
real(dp) :: t1, t2, t3
real(dp), allocatable :: x_real(:), wsave(:)
complex(dp), allocatable :: x(:), xdft(:), B(:)
call init_random()
n = 1024 * 1024
allocate(x_real(n), x(n), xdft(n), wsave(4*n+15), B(n*2))
print *, "FFTPACK"
call random_number(x_real)
x = x_real
call cpu_time(t1)
call zffti(n, wsave)
call cpu_time(t2)
call zfftf(n, x, wsave)
call cpu_time(t3)
!print *, x(:10)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "zffti:", (t2-t1)*1000, "ms"
print *, "zfftf:", (t3-t2)*1000, "ms"

print *, "NR"
x = x_real
call cpu_time(t1)
call dfour1(x, size(x), -1)
call cpu_time(t2)
!print *, x(:10)
print *, "Total time:", (t2-t1)*1000, "ms"

print *, "FFTE"
x = x_real
call cpu_time(t1)
call zfft1d(x, size(x), 0, B)
call cpu_time(t2)
call zfft1d(x, size(x), -1, B)
call cpu_time(t3)
!print *, x(:10)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "zfft1d init:", (t2-t1)*1000, "ms"
print *, "zfft1d:", (t3-t2)*1000, "ms"

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

end program
