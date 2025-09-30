      subroutine cosqb1(n,x,w,Xh)
      use fftpack_kind
      implicit none
      integer :: i , k , kc , modn , n , np2 , ns2
      real(rk) :: w , x , Xh , xim1
      dimension x(*) , w(*) , Xh(*)
      ns2 = (n+1)/2
      np2 = n + 2
      do i = 3 , n , 2
         xim1 = x(i-1) + x(i)
         x(i) = x(i) - x(i-1)
         x(i-1) = xim1
      enddo
      x(1) = x(1) + x(1)
      modn = mod(n,2)
      if ( modn==0 ) x(n) = x(n) + x(n)
      call dfftb(n,x,Xh)
      do k = 2 , ns2
         kc = np2 - k
         Xh(k) = w(k-1)*x(kc) + w(kc-1)*x(k)
         Xh(kc) = w(k-1)*x(k) - w(kc-1)*x(kc)
      enddo
      if ( modn==0 ) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))
      do k = 2 , ns2
         kc = np2 - k
         x(k) = Xh(k) + Xh(kc)
         x(kc) = Xh(k) - Xh(kc)
      enddo
      x(1) = x(1) + x(1)
      end subroutine cosqb1