program forward_transform_of_real_function
  !! This program invokes fftpack's rrft function to compute the forward transform of a real function
  !! and constructs the resulting complex Fourier coefficients by (re)organizing and normalizing the
  !! rfft result according to array element layout described at [1].  The program also demonstrates
  !! the inverse transform of the raw rrft result to recover the original function.
  !!
  !! [1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.rfft.html#scipy.fftpack.rfft
  use fftpack, only: rfft, irfft
  implicit none
  integer j, k
  integer, parameter :: N = 8
  double precision, parameter :: two_pi = 2.D0*acos(-1.D0), tolerance = 1.0D-06, f_avg = 3.D0, zero=0.D0
  double precision, parameter :: x(0:N-1) = [(two_pi*dble(j)/dble(N), j=0,N-1)]
  double precision, parameter :: f(0:N-1) = f_avg + cos(x)
    !! sample f(x) = 3 + cos(x) uniformly on [0,2*pi)
    !!             = 3 + (exp(i*x) - exp(-i*x))/2  
    !! which yields the Fourier coefficients
    !!                { 3, k = 0
    !!       f_hat =  { 1/2, k = 1
    !!                { 0, otherwise
  double precision, dimension(0:N-1) :: f_round_trip, rfft_f
  integer, parameter :: rk = kind(two_pi)
  complex(rk) f_hat(0:N/2)
  character(len=*), parameter :: real_format = "(a,*(g10.4,:,1x))" !! space-separated values
  character(len=*), parameter :: complex_format= "(a,*('(',g10.4,',',g10.4,')',:,1x)))" !! space-separated complex values

  call assert(mod(N,2)==0, "the algorithm below requires even N")

  rfft_f(:) = rfft(f)/dble(N)
  f_hat(:) = [ &
    cmplx(rfft_f(0),zero), &
    [( cmplx(rfft_f(k),rfft_f(k+1)),  k=lbound(rfft_f,1)+1,ubound(rfft_f,1)-1,2)], &
    cmplx(zero,rfft_f(N-1)) &
  ]
  f_round_trip(:) = irfft(rfft_f)

  print real_format, "f = ", f
  print complex_format, "f_hat = ", f_hat
  print real_format, "f_round_trip = ",f_round_trip

  call assert(any(abs(f_round_trip - f) < tolerance), "inverse of forward FFT must yield the original function")

contains

  pure subroutine assert(assertion, description)
    logical, intent(in) :: assertion
    character(len=*), intent(in) :: description
    if (.not. assertion) error stop description
  end subroutine

end program
