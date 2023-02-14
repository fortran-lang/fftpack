program main
use fftpack, only: zffti, zfftf, zfftb, fft, ifft
use fftpack_kind, only: rk

implicit none
integer, parameter :: N = 1024*135*77  !> (2**10)*(3**3)*5*7*11

complex(rk), dimension(N) :: x, z
real(rk), dimension(4*N+15) :: w
real(rk) :: err, time_i, time_f, time_b, t1, t2

call random_number(x%re)
z = x

print *, "02: Benchmarking zfft & fft"

call cpu_time(t1)
call zffti(N, w)
call cpu_time(t2)
time_i = t2-t1
print *, "Initializing: done"

call cpu_time(t1)
call zfftf(N, z, w)
call cpu_time(t2)
time_f = t2-t1
print *, "Forward: done"

call cpu_time(t1)
call zfftb(N, z, w)
call cpu_time(t2)
time_b = t2-t1
print *, "Backward: done"
print *, ""

err = maxval(abs(x-real(z/N,rk)))
print *, "--zfft"
print *, "Error: ", err
print *, "Init. time: ", time_i
print *, "Forward time: ", time_f
print *, "Backward time: ", time_b
print *, ""

print *, "Comparing to calls through fft"
call cpu_time(t1)
z = fft(x)
call cpu_time(t2)
time_f = t2-t1
print *, "Init. & forward: done"

call cpu_time(t1)
z = ifft(z)
call cpu_time(t2)
time_b = t2-t1
print *, "Backward: done"
print *, ""

err = maxval(abs(x-real(z/N,rk)))
print *, "--fft"
print *, "Error: ", err
print *, "Init. & forward time: ", time_f
print *, "Backward time: ", time_b
print *, ""
end program
