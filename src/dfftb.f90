      subroutine dfftb(n,r,Wsave)
      use fftpack_kind
      implicit none
      integer :: n
      real(rk) :: r , Wsave
      dimension r(*) , Wsave(*)
      if ( n==1 ) return
      call rfftb1(n,r,Wsave,Wsave(n+1),Wsave(2*n+1))
      end subroutine dfftb