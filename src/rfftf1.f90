      subroutine rfftf1(n, c, ch, wa, ifac)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(inout) :: c(*)
         real(dp), intent(in) :: wa(*)
         real(dp), intent(out) :: ch(*)
         integer, intent(in) :: ifac(*)
         integer :: i, idl1, ido, ip, iw, ix2, ix3, ix4, k1, &
                    kh, l1, l2, na, nf
         nf = ifac(2)
         na = 1
         l2 = n
         iw = n
         do k1 = 1, nf
            kh = nf - k1
            ip = ifac(kh + 3)
            l1 = l2/ip
            ido = n/l2
            idl1 = ido*l1
            iw = iw - (ip - 1)*ido
            na = 1 - na
            if (ip == 4) then
               ix2 = iw + ido
               ix3 = ix2 + ido
               if (na /= 0) then
                  call radf4(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
               else
                  call radf4(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
               end if
            elseif (ip /= 2) then
               if (ip == 3) then
                  ix2 = iw + ido
                  if (na /= 0) then
                     call radf3(ido, l1, ch, c, wa(iw), wa(ix2))
                  else
                     call radf3(ido, l1, c, ch, wa(iw), wa(ix2))
                  end if
               elseif (ip /= 5) then
                  if (ido == 1) na = 1 - na
                  if (na /= 0) then
                     call radfg(ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
                     na = 0
                  else
                     call radfg(ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
                     na = 1
                  end if
               else
                  ix2 = iw + ido
                  ix3 = ix2 + ido
                  ix4 = ix3 + ido
                  if (na /= 0) then
                     call radf5(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                  else
                     call radf5(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                  end if
               end if
            elseif (na /= 0) then
               call radf2(ido, l1, ch, c, wa(iw))
            else
               call radf2(ido, l1, c, ch, wa(iw))
            end if
            l2 = l1
         end do
         if (na == 1) return
         do i = 1, n
            c(i) = ch(i)
         end do
      end subroutine rfftf1
