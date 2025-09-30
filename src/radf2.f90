      subroutine radf2(ido, l1, cc, ch, wa1)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, l1, 2), wa1(*)
         real(dp), intent(out) :: ch(ido, 2, l1)
         real(dp) :: ti2, tr2
         integer :: i, ic, idp2, k
         do concurrent(k=1:l1)
            ch(1, 1, k) = cc(1, k, 1) + cc(1, k, 2)
            ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 2)
         end do
         if (ido < 2) return
         if (ido /= 2) then
            idp2 = ido + 2
            do concurrent(i=3:ido:2, k=1:l1)
               ic = idp2 - i
               tr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
               ti2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
               ch(i, 1, k) = cc(i, k, 1) + ti2
               ch(ic, 2, k) = ti2 - cc(i, k, 1)
               ch(i - 1, 1, k) = cc(i - 1, k, 1) + tr2
               ch(ic - 1, 2, k) = cc(i - 1, k, 1) - tr2
            end do
            if (mod(ido, 2) == 1) return
         end if
         do concurrent(k=1:l1)
            ch(1, 2, k) = -cc(ido, k, 2)
            ch(ido, 1, k) = cc(ido, k, 1)
         end do
      end subroutine radf2
