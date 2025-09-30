      subroutine cosqf1(n,x,w,Xh)
      use fftpack_kind
      implicit none
      integer :: i , k , kc , modn , n , np2 , ns2
      real(rk) :: w , x , Xh , xim1
      dimension x(*) , w(*) , Xh(*)
      ns2 = (n+1)/2
      np2 = n + 2
      do k = 2 , ns2
         kc = np2 - k
         Xh(k) = x(k) + x(kc)
         Xh(kc) = x(k) - x(kc)
      enddo
      modn = mod(n,2)
      if ( modn==0 ) Xh(ns2+1) = x(ns2+1) + x(ns2+1)
      do k = 2 , ns2
         kc = np2 - k
         x(k) = w(k-1)*Xh(kc) + w(kc-1)*Xh(k)
         x(kc) = w(k-1)*Xh(k) - w(kc-1)*Xh(kc)
      enddo
      if ( modn==0 ) x(ns2+1) = w(ns2)*Xh(ns2+1)
      call dfftf(n,x,Xh)
      do i = 3 , n , 2
         xim1 = x(i-1) - x(i)
         x(i) = x(i-1) + x(i)
         x(i-1) = xim1
      enddo
      end subroutine cosqf1