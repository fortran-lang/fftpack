      subroutine passf2(ido, l1, cc, ch, wa1)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 2, l1), wa1(*)
         real(dp), intent(out) :: ch(ido, l1, 2)
         real(dp) :: ti2, tr2
         integer :: i, k
         if (ido > 2) then
            do concurrent(k=1:l1, i=2:ido:2)
               ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(i - 1, 2, k)
               tr2 = cc(i - 1, 1, k) - cc(i - 1, 2, k)
               ch(i, k, 1) = cc(i, 1, k) + cc(i, 2, k)
               ti2 = cc(i, 1, k) - cc(i, 2, k)
               ch(i, k, 2) = wa1(i - 1)*ti2 - wa1(i)*tr2
               ch(i - 1, k, 2) = wa1(i - 1)*tr2 + wa1(i)*ti2
            end do
         else
            do concurrent(k=1:l1)
               ch(1, k, 1) = cc(1, 1, k) + cc(1, 2, k)
               ch(1, k, 2) = cc(1, 1, k) - cc(1, 2, k)
               ch(2, k, 1) = cc(2, 1, k) + cc(2, 2, k)
               ch(2, k, 2) = cc(2, 1, k) - cc(2, 2, k)
            end do
         end if
      end subroutine passf2
