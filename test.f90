program test
integer, parameter :: dp=kind(0.d0)
real(dp) :: t1, t2, t3
real(dp), allocatable :: x_real(:), wsave(:)
complex(dp), allocatable :: x(:), xdft(:)
print *, "test"
call init_random()
n = 1024 * 1024
allocate(x_real(n), x(n), xdft(n), wsave(4*n+15))
call random_number(x_real)
x = x_real
call cpu_time(t1)
call zffti(n, wsave)
call cpu_time(t2)
call zfftf(n, x, wsave)
call cpu_time(t3)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "zffti:", (t2-t1)*1000, "ms"
print *, "zfftf:", (t3-t2)*1000, "ms"

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

end program
