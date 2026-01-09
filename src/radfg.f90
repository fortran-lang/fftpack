      subroutine radfg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: ido, ip, l1, idl1
         real(dp), intent(in) :: wa(*)
         real(dp), intent(inout) :: cc(ido, ip, l1)
         real(dp), intent(inout) :: c1(ido, l1, ip)
         real(dp), intent(inout) :: c2(idl1, ip)
         real(dp), intent(out) :: ch(ido, l1, ip)
         real(dp), intent(inout) :: ch2(idl1, ip)
         real(dp) :: ai1, ai2, ar1, ar1h, ar2, ar2h, arg, &
                     dc2, dcp, ds2, dsp
         integer :: i, ic, idij, idp2, ik, ipp2, &
                    ipph, is, j, j2, jc, k, l, lc, nbd
         real(dp), parameter :: tpi = 2.0_dp*acos(-1.0_dp) ! 2 * pi
         arg = tpi/real(ip, kind=dp)
         dcp = cos(arg)
         dsp = sin(arg)
         ipph = (ip + 1)/2
         ipp2 = ip + 2
         idp2 = ido + 2
         nbd = (ido - 1)/2
         if (ido == 1) then
            do concurrent(ik=1:idl1)
               c2(ik, 1) = ch2(ik, 1)
            end do
         else
            do concurrent(ik=1:idl1)
               ch2(ik, 1) = c2(ik, 1)
            end do
            do concurrent(j=2:ip, k=1:l1)
               ch(1, k, j) = c1(1, k, j)
            end do
            if (nbd > l1) then
               is = -ido
               do j = 2, ip
                  is = is + ido
                  do k = 1, l1
                     idij = is
                     do i = 3, ido, 2
                        idij = idij + 2
                        ch(i - 1, k, j) = wa(idij - 1)*c1(i - 1, k, j) + wa(idij) &
                                          *c1(i, k, j)
                        ch(i, k, j) = wa(idij - 1)*c1(i, k, j) - wa(idij) &
                                      *c1(i - 1, k, j)
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
                        ch(i - 1, k, j) = wa(idij - 1)*c1(i - 1, k, j) + wa(idij) &
                                          *c1(i, k, j)
                        ch(i, k, j) = wa(idij - 1)*c1(i, k, j) - wa(idij) &
                                      *c1(i - 1, k, j)
                     end do
                  end do
               end do
            end if
            if (nbd < l1) then
               do concurrent(j=2:ipph, k=1:l1, i=3:ido:2)
                  jc = ipp2 - j
                  c1(i - 1, k, j) = ch(i - 1, k, j) + ch(i - 1, k, jc)
                  c1(i - 1, k, jc) = ch(i, k, j) - ch(i, k, jc)
                  c1(i, k, j) = ch(i, k, j) + ch(i, k, jc)
                  c1(i, k, jc) = ch(i - 1, k, jc) - ch(i - 1, k, j)
               end do
            else
               do concurrent(j=2:ipph, k=1:l1, i=3:ido:2)
                  jc = ipp2 - j
                  c1(i - 1, k, j) = ch(i - 1, k, j) + ch(i - 1, k, jc)
                  c1(i - 1, k, jc) = ch(i, k, j) - ch(i, k, jc)
                  c1(i, k, j) = ch(i, k, j) + ch(i, k, jc)
                  c1(i, k, jc) = ch(i - 1, k, jc) - ch(i - 1, k, j)
               end do
            end if
         end if
         do concurrent(j=2:ipph, k=1:l1)
            jc = ipp2 - j
            c1(1, k, j) = ch(1, k, j) + ch(1, k, jc)
            c1(1, k, jc) = ch(1, k, jc) - ch(1, k, j)
         end do
!
         ar1 = 1.0_dp
         ai1 = 0.0_dp
         do l = 2, ipph
            lc = ipp2 - l
            ar1h = dcp*ar1 - dsp*ai1
            ai1 = dcp*ai1 + dsp*ar1
            ar1 = ar1h
            do concurrent(ik=1:idl1)
               ch2(ik, l) = c2(ik, 1) + ar1*c2(ik, 2)
               ch2(ik, lc) = ai1*c2(ik, ip)
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
               do concurrent(ik=1:idl1)
                  ch2(ik, l) = ch2(ik, l) + ar2*c2(ik, j)
                  ch2(ik, lc) = ch2(ik, lc) + ai2*c2(ik, jc)
               end do
            end do
         end do
         do concurrent(j=2:ipph, ik=1:idl1)
            ch2(ik, 1) = ch2(ik, 1) + c2(ik, j)
         end do
!
         if (ido < l1) then
            do concurrent(k=1:l1, i=1:ido)
               cc(i, 1, k) = ch(i, k, 1)
            end do
         else
            do concurrent(i=1:ido, k=1:l1)
               cc(i, 1, k) = ch(i, k, 1)
            end do
         end if
         do concurrent(j=2:ipph, k=1:l1)
            jc = ipp2 - j
            j2 = j + j
            cc(ido, j2 - 2, k) = ch(1, k, j)
            cc(1, j2 - 1, k) = ch(1, k, jc)
         end do
         if (ido == 1) return
         if (nbd < l1) then
            do concurrent(j=2:ipph, k=1:l1, i=3:ido:2)
               jc = ipp2 - j
               j2 = j + j
               ic = idp2 - i
               cc(i - 1, j2 - 1, k) = ch(i - 1, k, j) + ch(i - 1, k, jc)
               cc(ic - 1, j2 - 2, k) = ch(i - 1, k, j) - ch(i - 1, k, jc)
               cc(i, j2 - 1, k) = ch(i, k, j) + ch(i, k, jc)
               cc(ic, j2 - 2, k) = ch(i, k, jc) - ch(i, k, j)
            end do
         else
            do concurrent(j=2:ipph, k=1:l1, i=3:ido:2)
               jc = ipp2 - j
               j2 = j + j
               ic = idp2 - i
               cc(i - 1, j2 - 1, k) = ch(i - 1, k, j) + ch(i - 1, k, jc)
               cc(ic - 1, j2 - 2, k) = ch(i - 1, k, j) - ch(i - 1, k, jc)
               cc(i, j2 - 1, k) = ch(i, k, j) + ch(i, k, jc)
               cc(ic, j2 - 2, k) = ch(i, k, jc) - ch(i, k, j)
            end do
         end if
      end subroutine radfg
