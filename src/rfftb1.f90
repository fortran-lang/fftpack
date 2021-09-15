!*==RFFTB1.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine rfftb1(n,c,Ch,Wa,Ifac)
      use fftpack_kind
      implicit none
!*--RFFTB12251
!*** Start of declarations inserted by SPAG
      real c , Ch , fftpack_kind , rk , Wa
      integer i , idl1 , ido , Ifac , ip , iw , ix2 , ix3 , ix4 , k1 ,  &
            & l1 , l2 , n , na , nf
!*** End of declarations inserted by SPAG
      dimension Ch(*) , c(*) , Wa(*) , Ifac(*)
      nf = Ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do k1 = 1 , nf
         ip = Ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if ( ip==4 ) then
            ix2 = iw + ido
            ix3 = ix2 + ido
            if ( na/=0 ) then
               call radb4(ido,l1,Ch,c,Wa(iw),Wa(ix2),Wa(ix3))
            else
               call radb4(ido,l1,c,Ch,Wa(iw),Wa(ix2),Wa(ix3))
            endif
            na = 1 - na
         elseif ( ip==2 ) then
            if ( na/=0 ) then
               call radb2(ido,l1,Ch,c,Wa(iw))
            else
               call radb2(ido,l1,c,Ch,Wa(iw))
            endif
            na = 1 - na
         elseif ( ip==3 ) then
            ix2 = iw + ido
            if ( na/=0 ) then
               call radb3(ido,l1,Ch,c,Wa(iw),Wa(ix2))
            else
               call radb3(ido,l1,c,Ch,Wa(iw),Wa(ix2))
            endif
            na = 1 - na
         elseif ( ip/=5 ) then
            if ( na/=0 ) then
               call radbg(ido,ip,l1,idl1,Ch,Ch,Ch,c,c,Wa(iw))
            else
               call radbg(ido,ip,l1,idl1,c,c,c,Ch,Ch,Wa(iw))
            endif
            if ( ido==1 ) na = 1 - na
         else
            ix2 = iw + ido
            ix3 = ix2 + ido
            ix4 = ix3 + ido
            if ( na/=0 ) then
               call radb5(ido,l1,Ch,c,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
            else
               call radb5(ido,l1,c,Ch,Wa(iw),Wa(ix2),Wa(ix3),Wa(ix4))
            endif
            na = 1 - na
         endif
         l1 = l2
         iw = iw + (ip-1)*ido
      enddo
      if ( na==0 ) return
      do i = 1 , n
         c(i) = Ch(i)
      enddo
      end subroutine rfftb1