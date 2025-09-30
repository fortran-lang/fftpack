      subroutine dffti(n,Wsave)
      use fftpack_kind
      implicit none
      integer :: n
      real(rk) :: Wsave
      dimension Wsave(*)
      if ( n==1 ) return
      call rffti1(n,Wsave(n+1),Wsave(2*n+1))
      end subroutine dffti