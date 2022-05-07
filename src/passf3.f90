      subroutine passf3(Ido,l1,Cc,Ch,Wa1,Wa2)
      use fftpack_kind
      implicit none
      real(rk) :: Cc , Ch , ci2 , ci3 , cr2 , cr3 , di2 , di3 , &
                     & dr2 , dr3 , ti2 , tr2 , Wa1 , Wa2
      integer :: i , Ido , k , l1
      dimension Cc(Ido,3,l1) , Ch(Ido,l1,3) , Wa1(*) , Wa2(*)
      real(rk),parameter :: taur = -0.5_rk
      real(rk),parameter :: taui = -sqrt(3.0_rk) / 2.0_rk
      if ( Ido/=2 ) then
         do k = 1 , l1
            do i = 2 , Ido , 2
               tr2 = Cc(i-1,2,k) + Cc(i-1,3,k)
               cr2 = Cc(i-1,1,k) + taur*tr2
               Ch(i-1,k,1) = Cc(i-1,1,k) + tr2
               ti2 = Cc(i,2,k) + Cc(i,3,k)
               ci2 = Cc(i,1,k) + taur*ti2
               Ch(i,k,1) = Cc(i,1,k) + ti2
               cr3 = taui*(Cc(i-1,2,k)-Cc(i-1,3,k))
               ci3 = taui*(Cc(i,2,k)-Cc(i,3,k))
               dr2 = cr2 - ci3
               dr3 = cr2 + ci3
               di2 = ci2 + cr3
               di3 = ci2 - cr3
               Ch(i,k,2) = Wa1(i-1)*di2 - Wa1(i)*dr2
               Ch(i-1,k,2) = Wa1(i-1)*dr2 + Wa1(i)*di2
               Ch(i,k,3) = Wa2(i-1)*di3 - Wa2(i)*dr3
               Ch(i-1,k,3) = Wa2(i-1)*dr3 + Wa2(i)*di3
            enddo
         enddo
      else
         do k = 1 , l1
            tr2 = Cc(1,2,k) + Cc(1,3,k)
            cr2 = Cc(1,1,k) + taur*tr2
            Ch(1,k,1) = Cc(1,1,k) + tr2
            ti2 = Cc(2,2,k) + Cc(2,3,k)
            ci2 = Cc(2,1,k) + taur*ti2
            Ch(2,k,1) = Cc(2,1,k) + ti2
            cr3 = taui*(Cc(1,2,k)-Cc(1,3,k))
            ci3 = taui*(Cc(2,2,k)-Cc(2,3,k))
            Ch(1,k,2) = cr2 - ci3
            Ch(1,k,3) = cr2 + ci3
            Ch(2,k,2) = ci2 + cr3
            Ch(2,k,3) = ci2 - cr3
         enddo
      end if
      end subroutine passf3