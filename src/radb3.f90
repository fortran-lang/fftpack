!*==RADB3.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine radb3(Ido,l1,Cc,Ch,Wa1,Wa2)
      use fftpack_kind
      implicit none
!*--RADB31487
!*** Start of declarations inserted by SPAG
      real Cc , Ch , ci2 , ci3 , cr2 , cr3 , di2 , di3 , dr2 , dr3 ,    &
         & fftpack_kind , rk , taui , taur , ti2 , tr2 , Wa1 , Wa2
      integer i , ic , Ido , idp2 , k , l1
!*** End of declarations inserted by SPAG
      dimension Cc(Ido,3,l1) , Ch(Ido,l1,3) , Wa1(1) , Wa2(1)
!     *** TAUI IS SQRT(3)/2 ***
      data taur , taui/ - 0.5d0 , 0.86602540378443864676d0/
      do k = 1 , l1
         tr2 = Cc(Ido,2,k) + Cc(Ido,2,k)
         cr2 = Cc(1,1,k) + taur*tr2
         Ch(1,k,1) = Cc(1,1,k) + tr2
         ci3 = taui*(Cc(1,3,k)+Cc(1,3,k))
         Ch(1,k,2) = cr2 - ci3
         Ch(1,k,3) = cr2 + ci3
      enddo
      if ( Ido==1 ) return
      idp2 = Ido + 2
      do k = 1 , l1
         do i = 3 , Ido , 2
            ic = idp2 - i
            tr2 = Cc(i-1,3,k) + Cc(ic-1,2,k)
            cr2 = Cc(i-1,1,k) + taur*tr2
            Ch(i-1,k,1) = Cc(i-1,1,k) + tr2
            ti2 = Cc(i,3,k) - Cc(ic,2,k)
            ci2 = Cc(i,1,k) + taur*ti2
            Ch(i,k,1) = Cc(i,1,k) + ti2
            cr3 = taui*(Cc(i-1,3,k)-Cc(ic-1,2,k))
            ci3 = taui*(Cc(i,3,k)+Cc(ic,2,k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            Ch(i-1,k,2) = Wa1(i-2)*dr2 - Wa1(i-1)*di2
            Ch(i,k,2) = Wa1(i-2)*di2 + Wa1(i-1)*dr2
            Ch(i-1,k,3) = Wa2(i-2)*dr3 - Wa2(i-1)*di3
            Ch(i,k,3) = Wa2(i-2)*di3 + Wa2(i-1)*dr3
         enddo
      enddo
      end subroutine radb3