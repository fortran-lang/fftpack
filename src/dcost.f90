      subroutine dcost(n,x,Wsave)
      use fftpack_kind
      implicit none
      real(rk) :: c1 , t1 , t2 , tx2 , Wsave , x , x1h , x1p3 , &
                  xi , xim2
      integer :: i , k , kc , modn , n , nm1 , np1 , ns2
      dimension x(*) , Wsave(*)
      nm1 = n - 1
      np1 = n + 1
      ns2 = n/2
      if ( n<2 ) return
      if ( n==2 ) then
         x1h = x(1) + x(2)
         x(2) = x(1) - x(2)
         x(1) = x1h
         return
      elseif ( n>3 ) then
         c1 = x(1) - x(n)
         x(1) = x(1) + x(n)
         do k = 2 , ns2
            kc = np1 - k
            t1 = x(k) + x(kc)
            t2 = x(k) - x(kc)
            c1 = c1 + Wsave(kc)*t2
            t2 = Wsave(k)*t2
            x(k) = t1 - t2
            x(kc) = t1 + t2
         enddo
         modn = mod(n,2)
         if ( modn/=0 ) x(ns2+1) = x(ns2+1) + x(ns2+1)
         call dfftf(nm1,x,Wsave(n+1))
         xim2 = x(2)
         x(2) = c1
         do i = 4 , n , 2
            xi = x(i)
            x(i) = x(i-2) - x(i-1)
            x(i-1) = xim2
            xim2 = xi
         enddo
         if ( modn/=0 ) x(n) = xim2
         return
      endif
      x1p3 = x(1) + x(3)
      tx2 = x(2) + x(2)
      x(2) = x(1) - x(3)
      x(1) = x1p3 + tx2
      x(3) = x1p3 - tx2
      end subroutine dcost