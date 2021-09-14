      subroutine zffti(n,Wsave)
      use fftpack_kind
      implicit none
      integer :: iw1 , iw2 , n
      real(rk) :: Wsave
      dimension Wsave(*)
      if ( n==1 ) return
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call cffti1(n,Wsave(iw1),Wsave(iw2))
      end subroutine zffti