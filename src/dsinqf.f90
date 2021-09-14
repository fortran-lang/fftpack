      subroutine dsinqf(n,x,Wsave)
      use fftpack_kind
      implicit none
      integer :: k , kc , n , ns2
      real(rk) :: Wsave , x , xhold
      dimension x(*) , Wsave(*)
      if ( n==1 ) return
      ns2 = n/2
      do k = 1 , ns2
         kc = n - k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
      enddo
      call dcosqf(n,x,Wsave)
      do k = 2 , n , 2
         x(k) = -x(k)
      enddo
      end subroutine dsinqf