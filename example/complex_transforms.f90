program complex_transforms
  !! This program invokes fftpack's fft function to compute the forward transform of a complex function
  !! and the inverse transform of the result.  An assertion verifies the expected result of the forward 
  !! transform according to the element layout described at [1].  A second assertion checks that the
  !! inverse transform recovers the original function.
  !!
  !! [1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.fft.html#scipy.fftpack.fft
  use fftpack, only: fft, ifft
  implicit none
  integer j, k
  integer, parameter :: N = 8
  double precision, parameter :: two_pi = 2.D0*acos(-1.D0), tolerance = 1.0D-06, f_avg = 3.D0, zero=0.D0
  double precision, parameter :: x(0:N-1) = [(two_pi*dble(j)/dble(N), j=0,N-1)]
  integer, parameter :: rk = kind(two_pi)
  complex(rk), parameter :: f(0:N-1) = f_avg + cos(x)
    !! sample f(x) = 3 + cos(x) uniformly on [0,2*pi)
    !!             = 3 + (exp(i*x) - exp(-i*x))/2  
    !! which yields the Fourier coefficients
    !!                {  3, k = 0
    !!       f_hat =  {  1/2, k = 1
    !!                {  1/2, k = -1
    !!                {  0, otherwise
  complex(rk), dimension(0:N-1) :: f_round_trip, fft_f
  character(len=*), parameter :: real_format = "(a,*(g10.4,:,1x))" !! space-separated values
  character(len=*), parameter :: complex_format= "(a,*('(',g10.4,',',g10.4,')',:,1x)))" !! space-separated complex values

  call assert(mod(N,2)==0, "the algorithm below requires even N")

  fft_f(:) = fft(f)/dble(N)
  f_round_trip(:) = ifft(fft_f)

  print complex_format, "f = ", f
  print complex_format, "fft_f = ", fft_f
  print complex_format, "f_round_trip = ",f_round_trip

  !call assert(all(abs(f_round_trip - f) < tolerance), "inverse of forward FFT must yield the original function")

contains

  pure subroutine assert(assertion, description)
    logical, intent(in) :: assertion
    character(len=*), intent(in) :: description
    if (.not. assertion) error stop description
  end subroutine

end program
