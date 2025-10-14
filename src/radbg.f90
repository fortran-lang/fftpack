      subroutine radbg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: ido, ip, l1, idl1
         real(dp), intent(in) :: cc(ido, ip, l1), wa(*)
         real(dp), intent(inout) :: c1(ido, l1, ip), ch2(idl1, ip)
         real(dp), intent(out) :: c2(idl1, ip), ch(ido, l1, ip)
         real(dp) :: ai1, ai2, ar1, ar1h, ar2, ar2h, arg, &
                     dc2, dcp, ds2, dsp
         integer :: i, ic, idij, idp2, ik, ipp2, &
                    ipph, is, j, j2, jc, k, l, lc, nbd
         real(dp), parameter :: tpi = 2*acos(-1.0_dp) ! 2 * pi
         arg = tpi/real(ip, dp)
         dcp = cos(arg)
         dsp = sin(arg)
         idp2 = ido + 2
         nbd = (ido - 1)/2
         ipp2 = ip + 2
         ipph = (ip + 1)/2
         if (ido < l1) then
            do i = 1, ido
               do k = 1, l1
                  ch(i, k, 1) = cc(i, 1, k)
               end do
            end do
         else
            do k = 1, l1
               do i = 1, ido
                  ch(i, k, 1) = cc(i, 1, k)
               end do
            end do
         end if
         do j = 2, ipph
            jc = ipp2 - j
            j2 = j + j
            do k = 1, l1
               ch(1, k, j) = cc(ido, j2 - 2, k) + cc(ido, j2 - 2, k)
               ch(1, k, jc) = cc(1, j2 - 1, k) + cc(1, j2 - 1, k)
            end do
         end do
         if (ido /= 1) then
            if (nbd < l1) then
               do j = 2, ipph
                  jc = ipp2 - j
                  do i = 3, ido, 2
                     ic = idp2 - i
                     do k = 1, l1
                        ch(i - 1, k, j) = cc(i - 1, 2*j - 1, k) + cc(ic - 1, 2*j - 2, k)
                        ch(i - 1, k, jc) = cc(i - 1, 2*j - 1, k) - cc(ic - 1, 2*j - 2, k)
                        ch(i, k, j) = cc(i, 2*j - 1, k) - cc(ic, 2*j - 2, k)
                        ch(i, k, jc) = cc(i, 2*j - 1, k) + cc(ic, 2*j - 2, k)
                     end do
                  end do
               end do
            else
               do j = 2, ipph
                  jc = ipp2 - j
                  do k = 1, l1
                     do i = 3, ido, 2
                        ic = idp2 - i
                        ch(i - 1, k, j) = cc(i - 1, 2*j - 1, k) + cc(ic - 1, 2*j - 2, k)
                        ch(i - 1, k, jc) = cc(i - 1, 2*j - 1, k) - cc(ic - 1, 2*j - 2, k)
                        ch(i, k, j) = cc(i, 2*j - 1, k) - cc(ic, 2*j - 2, k)
                        ch(i, k, jc) = cc(i, 2*j - 1, k) + cc(ic, 2*j - 2, k)
                     end do
                  end do
               end do
            end if
         end if
         ar1 = 1.0_dp
         ai1 = 0.0_dp
         do l = 2, ipph
            lc = ipp2 - l
            ar1h = dcp*ar1 - dsp*ai1
            ai1 = dcp*ai1 + dsp*ar1
            ar1 = ar1h
            do ik = 1, idl1
               c2(ik, l) = ch2(ik, 1) + ar1*ch2(ik, 2)
               c2(ik, lc) = ai1*ch2(ik, ip)
            end do
            dc2 = ar1
            ds2 = ai1
            ar2 = ar1
            ai2 = ai1
            do j = 3, ipph
               jc = ipp2 - j
               ar2h = dc2*ar2 - ds2*ai2
               ai2 = dc2*ai2 + ds2*ar2
               ar2 = ar2h
               do ik = 1, idl1
                  c2(ik, l) = c2(ik, l) + ar2*ch2(ik, j)
                  c2(ik, lc) = c2(ik, lc) + ai2*ch2(ik, jc)
               end do
            end do
         end do
         do j = 2, ipph
            do ik = 1, idl1
               ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
            end do
         end do
         do j = 2, ipph
            jc = ipp2 - j
            do k = 1, l1
               ch(1, k, j) = c1(1, k, j) - c1(1, k, jc)
               ch(1, k, jc) = c1(1, k, j) + c1(1, k, jc)
            end do
         end do
         if (ido /= 1) then
            if (nbd < l1) then
               do j = 2, ipph
                  jc = ipp2 - j
                  do i = 3, ido, 2
                     do k = 1, l1
                        ch(i - 1, k, j) = c1(i - 1, k, j) - c1(i, k, jc)
                        ch(i - 1, k, jc) = c1(i - 1, k, j) + c1(i, k, jc)
                        ch(i, k, j) = c1(i, k, j) + c1(i - 1, k, jc)
                        ch(i, k, jc) = c1(i, k, j) - c1(i - 1, k, jc)
                     end do
                  end do
               end do
            else
               do j = 2, ipph
                  jc = ipp2 - j
                  do k = 1, l1
                     do i = 3, ido, 2
                        ch(i - 1, k, j) = c1(i - 1, k, j) - c1(i, k, jc)
                        ch(i - 1, k, jc) = c1(i - 1, k, j) + c1(i, k, jc)
                        ch(i, k, j) = c1(i, k, j) + c1(i - 1, k, jc)
                        ch(i, k, jc) = c1(i, k, j) - c1(i - 1, k, jc)
                     end do
                  end do
               end do
            end if
         end if
         if (ido == 1) return
         do ik = 1, idl1
            c2(ik, 1) = ch2(ik, 1)
         end do
         do j = 2, ip
            do k = 1, l1
               c1(1, k, j) = ch(1, k, j)
            end do
         end do
         if (nbd > l1) then
            is = -ido
            do j = 2, ip
               is = is + ido
               do k = 1, l1
                  idij = is
                  do i = 3, ido, 2
                     idij = idij + 2
                     c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij) &
                                       *ch(i, k, j)
                     c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij) &
                                   *ch(i - 1, k, j)
                  end do
               end do
            end do
         else
            is = -ido
            do j = 2, ip
               is = is + ido
               idij = is
               do i = 3, ido, 2
                  idij = idij + 2
                  do k = 1, l1
                     c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij) &
                                       *ch(i, k, j)
                     c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij) &
                                   *ch(i - 1, k, j)
                  end do
               end do
            end do
         end if
      end subroutine radbg
