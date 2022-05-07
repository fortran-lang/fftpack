      subroutine dsint(n,x,Wsave)
      use fftpack_kind
      implicit none
      integer :: iw1 , iw2 , iw3 , n , np1
      real(rk) :: Wsave , x
      dimension x(*) , Wsave(*)
      np1 = n + 1
      iw1 = n/2 + 1
      iw2 = iw1 + np1
      iw3 = iw2 + np1
      call sint1(n,x,Wsave,Wsave(iw1),Wsave(iw2),Wsave(iw3))
      end subroutine dsint