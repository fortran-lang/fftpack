program forward_transform_of_real_function
  !! This program computes the forward transform of a real function and constructs
  !! the complex result (re)organized to match the array subscripting to the common
  !! analytical form [1].  Which form one uses in practice requires balancing the
  !! need for speed versus clarity. 
  !!
  !! [1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.rfft.html#scipy.fftpack.rfft
  use fftpack, only: rfft, irfft
  implicit none
  integer j, k
  integer, parameter :: N = 8
  double precision, parameter :: two_pi = 2.D0*acos(-1.D0), tolerance = 1.0D-06, f_avg = 3.D0, zero=0.D0
  double precision, parameter :: f(0:N-1) = f_avg + [(cos(two_pi*dble(j)/dble(N)), j=0,N-1)]
  double precision, dimension(0:N-1) :: f_round_trip, rfft_f
  integer, parameter :: rk = kind(two_pi)
  complex(rk) f_hat(0:N/2)

  call assert(mod(N,2)==0, "the algorithm below requires even N")

  rfft_f(:) = rfft(f)/dble(N)
  f_hat(:) = [ cmplx(rfft_f(0),zero), [( cmplx(k,k+1),  k=lbound(rfft_f,1)+1,ubound(rfft_f,1)-1,2)], cmplx(zero,rfft_f(N-1)) ]
  f_round_trip(:) = dble(irfft(rfft_f))
  !call assert(any(abs(f_round_trip - f) < tolerance), "inverse of forward FFT must yield the original function")

  print *, "f = ", f
  print *, "f_hat = ", f_hat
  print *, "f_round_trip = ", f_round_trip

  !print '(3(10x,a,10x))',"f", "f_round_trip", "rfft_f"
  !do m = 1, size(f)
  !  print *, f(m), f_round_trip(m), rfft_f(m)
  !end do
  !print *

contains
  pure subroutine assert(assertion, description)
    logical, intent(in) :: assertion
    character(len=*), intent(in) :: description
    if (.not. assertion) error stop description
  end subroutine
end program
