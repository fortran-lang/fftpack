program bench1
use fftpack, only: zffti, zfftf, zfftb
implicit none
integer, parameter :: dp = kind(0.d0)
complex(dp), allocatable :: z(:)
real(dp), allocatable :: w(:), x(:)
real(dp) :: err, time_init, time_forward, time_backward, t1, t2
integer :: N

N = 1024*1014*16

allocate(x(N), z(N), w(4*N+15))
call random_number(x)
z = x

print *, "Initializing"
call cpu_time(t1)
call zffti(N, w)
call cpu_time(t2)
time_init = t2-t1

print *, "Forward"
call cpu_time(t1)
call zfftf(N, z, w)
call cpu_time(t2)
time_forward = t2-t1

print *, "Backward"
call cpu_time(t1)
call zfftb(N, z, w)
call cpu_time(t2)
time_backward = t2-t1
print *, "Done"

err = maxval(abs(x-real(z/N,dp)))
print *
print *, "Error: ", err
print *, "Init time: ", time_init
print *, "Forward time: ", time_forward
print *, "Backward time: ", time_backward
end program
