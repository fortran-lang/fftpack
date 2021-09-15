!*==RFFTF1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine rfftf1(n,c,Ch,Wa,Ifac)
      use fftpack_kind
      implicit none
!*--RFFTF12321
!*** Start of declarations inserted by SPAG
      real c , Ch , fftpack_kind , rk , Wa
      integer i , idl1 , ido , Ifac , ip , iw , ix2 , ix3 , ix4 , k1 ,  &
            & kh , l1 , l2 , n , na , nf
!*** End of declarations inserted by SPAG
      dimension Ch(*) , c(*) , Wa(*) , Ifac(*)
      nf = Ifac(2)
      na = 1
      l2 = n
      iw = n
      do k1 = 1 , nf
         kh = nf - k1
         ip = Ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw - (ip-1)*ido
         na = 1 - na
         if ( ip==4 ) then
            ix2 = iw + ido
            ix3 = ix2 + ido
            if ( na/=0 ) then
               call radf4(ido,l1,Ch,c,Wa(iw),Wa(ix2),Wa(ix3))
            else
               call radf4(ido,l1,c,Ch,Wa(iw),Wa(ix2),Wa(ix3))
            endif
         elseif ( ip/=2 ) then
            if ( ip==3 ) then
               ix2 = iw + ido
               if ( na/=0 ) then
                  call radf3(ido,l1,Ch,c,Wa(iw),Wa(ix2))
               else
                  call radf3(ido,l1,c,Ch,Wa(iw),Wa(ix2))
               endif
            elseif ( ip/=5 ) then
               if ( ido==1 ) na = 1 - na
               if ( na/=0 ) then
                  call radfg(ido,ip,l1,idl1,Ch,Ch,Ch,c,c,Wa(iw))
                  na = 0
               else
                  call radfg(ido,ip,l1,idl1,c,c,c,Ch,Ch,Wa(iw))
                  na = 1
               endif
            else
               ix2 = iw + ido
               ix3 = ix2 + ido
               ix4 = ix3 + ido
               if ( na/=0 ) then
                  call radf5(ido,l1,Ch,c,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
               else
                  call radf5(ido,l1,c,Ch,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
               endif
            endif
         elseif ( na/=0 ) then
            call radf2(ido,l1,Ch,c,Wa(iw))
         else
            call radf2(ido,l1,c,Ch,Wa(iw))
         endif
         l2 = l1
      enddo
      if ( na==1 ) return
      do i = 1 , n
         c(i) = Ch(i)
      enddo
      end subroutine rfftf1