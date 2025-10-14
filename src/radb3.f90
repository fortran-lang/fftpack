      subroutine radb3(ido, l1, cc, ch, wa1, wa2)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 3, l1), wa1(*), wa2(*)
         real(dp), intent(out) :: ch(ido, l1, 3)
         real(dp) :: ci2, ci3, cr2, cr3, di2, di3, &
                     dr2, dr3, ti2, tr2
         integer :: i, ic, idp2, k
         real(dp), parameter :: taur = -0.5_dp
         real(dp), parameter :: taui = sqrt(3.0_dp)/2.0_dp
         do k = 1, l1
            tr2 = cc(ido, 2, k) + cc(ido, 2, k)
            cr2 = cc(1, 1, k) + taur*tr2
            ch(1, k, 1) = cc(1, 1, k) + tr2
            ci3 = taui*(cc(1, 3, k) + cc(1, 3, k))
            ch(1, k, 2) = cr2 - ci3
            ch(1, k, 3) = cr2 + ci3
         end do
         if (ido == 1) return
         idp2 = ido + 2
         do k = 1, l1
            do i = 3, ido, 2
               ic = idp2 - i
               tr2 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
               cr2 = cc(i - 1, 1, k) + taur*tr2
               ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
               ti2 = cc(i, 3, k) - cc(ic, 2, k)
               ci2 = cc(i, 1, k) + taur*ti2
               ch(i, k, 1) = cc(i, 1, k) + ti2
               cr3 = taui*(cc(i - 1, 3, k) - cc(ic - 1, 2, k))
               ci3 = taui*(cc(i, 3, k) + cc(ic, 2, k))
               dr2 = cr2 - ci3
               dr3 = cr2 + ci3
               di2 = ci2 + cr3
               di3 = ci2 - cr3
               ch(i - 1, k, 2) = wa1(i - 2)*dr2 - wa1(i - 1)*di2
               ch(i, k, 2) = wa1(i - 2)*di2 + wa1(i - 1)*dr2
               ch(i - 1, k, 3) = wa2(i - 2)*dr3 - wa2(i - 1)*di3
               ch(i, k, 3) = wa2(i - 2)*di3 + wa2(i - 1)*dr3
            end do
         end do
      end subroutine radb3
