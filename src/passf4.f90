      subroutine passf4(Ido,l1,Cc,Ch,Wa1,Wa2,Wa3)
      use fftpack_kind
      implicit none
      real(rk) :: Cc , Ch , ci2 , ci3 , ci4 , cr2 , cr3 , cr4 , &
                     & ti1 , ti2 , ti3 , ti4 , tr1 , tr2 , tr3 , tr4 , &
                     & Wa1 , Wa2 , Wa3
      integer :: i , Ido , k , l1
      dimension Cc(Ido,4,l1) , Ch(Ido,l1,4) , Wa1(*) , Wa2(*) , Wa3(*)
      if ( Ido/=2 ) then
         do k = 1 , l1
            do i = 2 , Ido , 2
               ti1 = Cc(i,1,k) - Cc(i,3,k)
               ti2 = Cc(i,1,k) + Cc(i,3,k)
               ti3 = Cc(i,2,k) + Cc(i,4,k)
               tr4 = Cc(i,2,k) - Cc(i,4,k)
               tr1 = Cc(i-1,1,k) - Cc(i-1,3,k)
               tr2 = Cc(i-1,1,k) + Cc(i-1,3,k)
               ti4 = Cc(i-1,4,k) - Cc(i-1,2,k)
               tr3 = Cc(i-1,2,k) + Cc(i-1,4,k)
               Ch(i-1,k,1) = tr2 + tr3
               cr3 = tr2 - tr3
               Ch(i,k,1) = ti2 + ti3
               ci3 = ti2 - ti3
               cr2 = tr1 + tr4
               cr4 = tr1 - tr4
               ci2 = ti1 + ti4
               ci4 = ti1 - ti4
               Ch(i-1,k,2) = Wa1(i-1)*cr2 + Wa1(i)*ci2
               Ch(i,k,2) = Wa1(i-1)*ci2 - Wa1(i)*cr2
               Ch(i-1,k,3) = Wa2(i-1)*cr3 + Wa2(i)*ci3
               Ch(i,k,3) = Wa2(i-1)*ci3 - Wa2(i)*cr3
               Ch(i-1,k,4) = Wa3(i-1)*cr4 + Wa3(i)*ci4
               Ch(i,k,4) = Wa3(i-1)*ci4 - Wa3(i)*cr4
            enddo
         enddo
      else
         do k = 1 , l1
            ti1 = Cc(2,1,k) - Cc(2,3,k)
            ti2 = Cc(2,1,k) + Cc(2,3,k)
            tr4 = Cc(2,2,k) - Cc(2,4,k)
            ti3 = Cc(2,2,k) + Cc(2,4,k)
            tr1 = Cc(1,1,k) - Cc(1,3,k)
            tr2 = Cc(1,1,k) + Cc(1,3,k)
            ti4 = Cc(1,4,k) - Cc(1,2,k)
            tr3 = Cc(1,2,k) + Cc(1,4,k)
            Ch(1,k,1) = tr2 + tr3
            Ch(1,k,3) = tr2 - tr3
            Ch(2,k,1) = ti2 + ti3
            Ch(2,k,3) = ti2 - ti3
            Ch(1,k,2) = tr1 + tr4
            Ch(1,k,4) = tr1 - tr4
            Ch(2,k,2) = ti1 + ti4
            Ch(2,k,4) = ti1 - ti4
         enddo
      end if
      end subroutine passf4