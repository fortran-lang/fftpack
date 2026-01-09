      subroutine passf(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(out) :: nac
         integer, intent(in) :: ido, ip, l1, idl1
         real(dp), intent(in) :: cc(ido, ip, l1), wa(*)
         real(dp), intent(out) :: c1(ido, l1, ip), c2(idl1, ip), ch(ido, l1, ip)
         real(dp), intent(inout) :: ch2(idl1, ip)
         real(dp) :: wai, war
         integer :: i, idij, idj, idl, idlj, idot, idp, &
                    ik, inc, ipp2, ipph, j, jc, k, l, lc
         integer :: nt
         idot = ido/2
         nt = ip*idl1
         ipp2 = ip + 2
         ipph = (ip + 1)/2
         idp = ip*ido
!
         if (ido < l1) then
            do concurrent(k=1:l1, j=2:ipph, i=1:ido)
               jc = ipp2 - j
               ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
               ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
            end do
            do concurrent(k=1:l1, i=1:ido)
               ch(i, k, 1) = cc(i, 1, k)
            end do
         else
            do concurrent(k=1:l1, j=2:ipph, i=1:ido)
               jc = ipp2 - j
               ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
               ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
            end do
            do concurrent(i=1:ido, k=1:l1)
               ch(i, k, 1) = cc(i, 1, k)
            end do
         end if
         idl = 2 - ido
         inc = 0
         do l = 2, ipph
            lc = ipp2 - l
            idl = idl + ido
            do concurrent(ik=1:idl1)
               c2(ik, l) = ch2(ik, 1) + wa(idl - 1)*ch2(ik, 2)
               c2(ik, lc) = -wa(idl)*ch2(ik, ip)
            end do
            idlj = idl
            inc = inc + ido
            do j = 3, ipph
               jc = ipp2 - j
               idlj = idlj + inc
               if (idlj > idp) idlj = idlj - idp
               war = wa(idlj - 1)
               wai = wa(idlj)
               do concurrent(ik=1:idl1)
                  c2(ik, l) = c2(ik, l) + war*ch2(ik, j)
                  c2(ik, lc) = c2(ik, lc) - wai*ch2(ik, jc)
               end do
            end do
         end do
         do concurrent(ik=1:idl1, j=2:ipph)
            ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
         end do
         do concurrent(j=2:ipph, ik=2:idl1:2)
            jc = ipp2 - j
            ch2(ik - 1, j) = c2(ik - 1, j) - c2(ik, jc)
            ch2(ik - 1, jc) = c2(ik - 1, j) + c2(ik, jc)
            ch2(ik, j) = c2(ik, j) + c2(ik - 1, jc)
            ch2(ik, jc) = c2(ik, j) - c2(ik - 1, jc)
         end do
         nac = 1
         if (ido == 2) return
         nac = 0
         do concurrent(ik=1:idl1)
            c2(ik, 1) = ch2(ik, 1)
         end do
         do concurrent(j=2:ip, k=1:l1)
            c1(1, k, j) = ch(1, k, j)
            c1(2, k, j) = ch(2, k, j)
         end do
         if (idot > l1) then
            idj = 2 - ido
            do j = 2, ip
               idj = idj + ido
               do k = 1, l1
                  idij = idj
                  do i = 4, ido, 2
                     idij = idij + 2
                     c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) + wa(idij) &
                                       *ch(i, k, j)
                     c1(i, k, j) = wa(idij - 1)*ch(i, k, j) - wa(idij) &
                                   *ch(i - 1, k, j)
                  end do
               end do
            end do
         else
            idij = 0
            do j = 2, ip
               idij = idij + 2
               do i = 4, ido, 2
                  idij = idij + 2
                  do concurrent(k=1:l1)
                     c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) + wa(idij)*ch(i, k, j)
                     c1(i, k, j) = wa(idij - 1)*ch(i, k, j) - wa(idij)*ch(i - 1, k, j)
                  end do
               end do
            end do
         end if
      end subroutine passf
