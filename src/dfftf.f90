      subroutine dfftf(n,r,Wsave)
      use fftpack_kind
      implicit none
      integer :: n
      real(rk) :: r , Wsave
      dimension r(1) , Wsave(*)
      if ( n==1 ) return
      call rfftf1(n,r,Wsave,Wsave(n+1),Wsave(2*n+1))
      end subroutine dfftf