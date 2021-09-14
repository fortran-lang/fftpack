      subroutine dsinqb(n,x,Wsave)
      use fftpack_kind
      implicit none
      integer :: k , kc , n , ns2
      real(rk) :: Wsave , x , xhold
      dimension x(*) , Wsave(*)
      if ( n>1 ) then
         ns2 = n/2
         do k = 2 , n , 2
            x(k) = -x(k)
         enddo
         call dcosqb(n,x,Wsave)
         do k = 1 , ns2
            kc = n - k
            xhold = x(k)
            x(k) = x(kc+1)
            x(kc+1) = xhold
         enddo
         return
      endif
      x(1) = 4.0_rk*x(1)
      return
      end subroutine dsinqb