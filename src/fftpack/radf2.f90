      subroutine radf2(Ido,l1,Cc,Ch,Wa1)
      use fftpack_kind
      implicit none
      real(rk) :: Cc , Ch , ti2 , tr2 , Wa1
      integer :: i , ic , Ido , idp2 , k , l1
      dimension Ch(Ido,2,l1) , Cc(Ido,l1,2) , Wa1(*)
      do k = 1 , l1
         Ch(1,1,k) = Cc(1,k,1) + Cc(1,k,2)
         Ch(Ido,2,k) = Cc(1,k,1) - Cc(1,k,2)
      enddo
      if ( Ido<2 ) return
      if ( Ido/=2 ) then
         idp2 = Ido + 2
         do k = 1 , l1
            do i = 3 , Ido , 2
               ic = idp2 - i
               tr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
               ti2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
               Ch(i,1,k) = Cc(i,k,1) + ti2
               Ch(ic,2,k) = ti2 - Cc(i,k,1)
               Ch(i-1,1,k) = Cc(i-1,k,1) + tr2
               Ch(ic-1,2,k) = Cc(i-1,k,1) - tr2
            enddo
         enddo
         if ( mod(Ido,2)==1 ) return
      endif
      do k = 1 , l1
         Ch(1,2,k) = -Cc(Ido,k,2)
         Ch(Ido,1,k) = Cc(Ido,k,1)
      enddo
      end subroutine radf2