      subroutine sint1(n,War,Was,Xh,x,Ifac)
      use fftpack_kind
      implicit none
      integer :: i , Ifac , k , kc , modn , n , np1 , ns2
      real(rk) :: t1 , t2 , War , Was , x , Xh , xhold
      dimension War(*) , Was(*) , x(*) , Xh(*) , Ifac(*)
      real(rk),parameter :: sqrt3 = sqrt(3.0_rk)
      do i = 1 , n
         Xh(i) = War(i)
         War(i) = x(i)
      enddo
      if ( n<2 ) then
         Xh(1) = Xh(1) + Xh(1)
      elseif ( n==2 ) then
         xhold = sqrt3*(Xh(1)+Xh(2))
         Xh(2) = sqrt3*(Xh(1)-Xh(2))
         Xh(1) = xhold
      else
         np1 = n + 1
         ns2 = n/2
         x(1) = 0.0_rk
         do k = 1 , ns2
            kc = np1 - k
            t1 = Xh(k) - Xh(kc)
            t2 = Was(k)*(Xh(k)+Xh(kc))
            x(k+1) = t1 + t2
            x(kc+1) = t2 - t1
         enddo
         modn = mod(n,2)
         if ( modn/=0 ) x(ns2+2) = 4.0_rk*Xh(ns2+1)
         call rfftf1(np1,x,Xh,War,Ifac)
         Xh(1) = 0.5_rk*x(1)
         do i = 3 , n , 2
            Xh(i-1) = -x(i)
            Xh(i) = Xh(i-2) + x(i-1)
         enddo
         if ( modn==0 ) Xh(n) = -x(n+1)
      endif
      do i = 1 , n
         x(i) = War(i)
         War(i) = Xh(i)
      enddo
      end subroutine sint1