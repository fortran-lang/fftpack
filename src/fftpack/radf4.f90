      subroutine radf4(Ido,l1,Cc,Ch,Wa1,Wa2,Wa3)
      use fftpack_kind
      implicit none
      real(rk) :: Cc , Ch , ci2 , ci3 , ci4 , cr2 , cr3 , cr4 , &
                  ti1 , ti2 , ti3 , ti4 , tr1 , tr2 , tr3, &
                  tr4 , Wa1 , Wa2 , Wa3
      integer :: i , ic , Ido , idp2 , k , l1
      dimension Cc(Ido,l1,4) , Ch(Ido,4,l1) , Wa1(*) , Wa2(*) , Wa3(*)
      real(rk),parameter :: hsqt2 = sqrt(2.0_rk) / 2.0_rk
      do k = 1 , l1
         tr1 = Cc(1,k,2) + Cc(1,k,4)
         tr2 = Cc(1,k,1) + Cc(1,k,3)
         Ch(1,1,k) = tr1 + tr2
         Ch(Ido,4,k) = tr2 - tr1
         Ch(Ido,2,k) = Cc(1,k,1) - Cc(1,k,3)
         Ch(1,3,k) = Cc(1,k,4) - Cc(1,k,2)
      enddo
      if ( Ido<2 ) return
      if ( Ido/=2 ) then
         idp2 = Ido + 2
         do k = 1 , l1
            do i = 3 , Ido , 2
               ic = idp2 - i
               cr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
               ci2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
               cr3 = Wa2(i-2)*Cc(i-1,k,3) + Wa2(i-1)*Cc(i,k,3)
               ci3 = Wa2(i-2)*Cc(i,k,3) - Wa2(i-1)*Cc(i-1,k,3)
               cr4 = Wa3(i-2)*Cc(i-1,k,4) + Wa3(i-1)*Cc(i,k,4)
               ci4 = Wa3(i-2)*Cc(i,k,4) - Wa3(i-1)*Cc(i-1,k,4)
               tr1 = cr2 + cr4
               tr4 = cr4 - cr2
               ti1 = ci2 + ci4
               ti4 = ci2 - ci4
               ti2 = Cc(i,k,1) + ci3
               ti3 = Cc(i,k,1) - ci3
               tr2 = Cc(i-1,k,1) + cr3
               tr3 = Cc(i-1,k,1) - cr3
               Ch(i-1,1,k) = tr1 + tr2
               Ch(ic-1,4,k) = tr2 - tr1
               Ch(i,1,k) = ti1 + ti2
               Ch(ic,4,k) = ti1 - ti2
               Ch(i-1,3,k) = ti4 + tr3
               Ch(ic-1,2,k) = tr3 - ti4
               Ch(i,3,k) = tr4 + ti3
               Ch(ic,2,k) = tr4 - ti3
            enddo
         enddo
         if ( mod(Ido,2)==1 ) return
      endif
      do k = 1 , l1
         ti1 = -hsqt2*(Cc(Ido,k,2)+Cc(Ido,k,4))
         tr1 = hsqt2*(Cc(Ido,k,2)-Cc(Ido,k,4))
         Ch(Ido,1,k) = tr1 + Cc(Ido,k,1)
         Ch(Ido,3,k) = Cc(Ido,k,1) - tr1
         Ch(1,2,k) = ti1 - Cc(Ido,k,3)
         Ch(1,4,k) = ti1 + Cc(Ido,k,3)
      enddo
      end subroutine radf4