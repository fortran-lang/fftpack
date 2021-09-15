!*==RADF3.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine radf3(Ido,l1,Cc,Ch,Wa1,Wa2)
      use fftpack_kind
      implicit none
!*--RADF31889
!*** Start of declarations inserted by SPAG
      real Cc , Ch , ci2 , cr2 , di2 , di3 , dr2 , dr3 , fftpack_kind , &
         & rk , taui , taur , ti2 , ti3 , tr2 , tr3 , Wa1 , Wa2
      integer i , ic , Ido , idp2 , k , l1
!*** End of declarations inserted by SPAG
      dimension Ch(Ido,3,l1) , Cc(Ido,l1,3) , Wa1(1) , Wa2(1)
!     *** TAUI IS -SQRT(3)/2 ***
      data taur , taui/ - 0.5d0 , 0.86602540378443864676d0/
      do k = 1 , l1
         cr2 = Cc(1,k,2) + Cc(1,k,3)
         Ch(1,1,k) = Cc(1,k,1) + cr2
         Ch(1,3,k) = taui*(Cc(1,k,3)-Cc(1,k,2))
         Ch(Ido,2,k) = Cc(1,k,1) + taur*cr2
      enddo
      if ( Ido==1 ) return
      idp2 = Ido + 2
      do k = 1 , l1
         do i = 3 , Ido , 2
            ic = idp2 - i
            dr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
            di2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
            dr3 = Wa2(i-2)*Cc(i-1,k,3) + Wa2(i-1)*Cc(i,k,3)
            di3 = Wa2(i-2)*Cc(i,k,3) - Wa2(i-1)*Cc(i-1,k,3)
            cr2 = dr2 + dr3
            ci2 = di2 + di3
            Ch(i-1,1,k) = Cc(i-1,k,1) + cr2
            Ch(i,1,k) = Cc(i,k,1) + ci2
            tr2 = Cc(i-1,k,1) + taur*cr2
            ti2 = Cc(i,k,1) + taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            Ch(i-1,3,k) = tr2 + tr3
            Ch(ic-1,2,k) = tr2 - tr3
            Ch(i,3,k) = ti2 + ti3
            Ch(ic,2,k) = ti3 - ti2
         enddo
      enddo
      end subroutine radf3