program bench3
   use fftpack, only: dffti, dfftf, dfftb, rfft, irfft
   use fftpack_kinds, only: dp

   implicit none
   integer, parameter :: N = 1024*135*77  ! (2**10)*(3**3)*5*7*11

   real(dp), dimension(N) :: x, y
   real(dp), dimension(2*N + 15) :: w
   real(dp) :: err, time_i, time_f, time_b, t1, t2

   call random_number(x)
   y = x

   print *, "03: Benchmarking dfft & rfft"

   call cpu_time(t1)
   call dffti(N, w)
   call cpu_time(t2)
   time_i = t2 - t1
   print *, "Initializing: done"

   call cpu_time(t1)
   call dfftf(N, y, w)
   call cpu_time(t2)
   time_f = t2 - t1
   print *, "Forward: done"

   call cpu_time(t1)
   call dfftb(N, y, w)
   call cpu_time(t2)
   time_b = t2 - t1
   print *, "Backward: done"
   print *, ""

   err = maxval(abs(x - y/N))
   print *, "--dfft"
   print *, "Error: ", err
   print *, "Init. time: ", time_i
   print *, "Forward time: ", time_f
   print *, "Backward time: ", time_b
   print *, ""

   print *, "Comparing to calls through rfft"
   call cpu_time(t1)
   y = rfft(x)
   call cpu_time(t2)
   time_f = t2 - t1
   print *, "Init. & forward: done"

   call cpu_time(t1)
   y = irfft(y)
   call cpu_time(t2)
   time_b = t2 - t1
   print *, "Backward: done"
   print *, ""

   err = maxval(abs(x - y/N))
   print *, "--rfft"
   print *, "Error: ", err
   print *, "Init. & forward time: ", time_f
   print *, "Backward time: ", time_b
   print *, ""
end program
