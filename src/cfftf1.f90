      subroutine cfftf1(n, c, ch, wa, ifac)
         use fftpack_kinds, only: dp
         use fftpack_legacy_drivers_pass, only: passf, passf2, passf3, passf4, passf5
         implicit none
         integer, intent(in) :: n, ifac(*)
         real(dp), intent(inout) :: c(*), ch(*)
         real(dp), intent(in) :: wa(*)
         integer :: i, idl1, ido, idot, ip, iw, ix2, ix3, ix4, &
                    k1, l1, l2, n2, na, nac, nf
         nf = ifac(2)
         na = 0
         l1 = 1
         iw = 1
         do k1 = 1, nf
            ip = ifac(k1 + 2)
            l2 = ip*l1
            ido = n/l2
            idot = ido + ido
            idl1 = idot*l1
            if (ip == 4) then
               ix2 = iw + idot
               ix3 = ix2 + idot
               if (na /= 0) then
                  call passf4(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
               else
                  call passf4(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
               end if
               na = 1 - na
            elseif (ip == 2) then
               if (na /= 0) then
                  call passf2(idot, l1, ch, c, wa(iw))
               else
                  call passf2(idot, l1, c, ch, wa(iw))
               end if
               na = 1 - na
            elseif (ip == 3) then
               ix2 = iw + idot
               if (na /= 0) then
                  call passf3(idot, l1, ch, c, wa(iw), wa(ix2))
               else
                  call passf3(idot, l1, c, ch, wa(iw), wa(ix2))
               end if
               na = 1 - na
            elseif (ip /= 5) then
               if (na /= 0) then
                  call passf(nac, idot, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
               else
                  call passf(nac, idot, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
               end if
               if (nac /= 0) na = 1 - na
            else
               ix2 = iw + idot
               ix3 = ix2 + idot
               ix4 = ix3 + idot
               if (na /= 0) then
                  call passf5(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
               else
                  call passf5(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
               end if
               na = 1 - na
            end if
            l1 = l2
            iw = iw + (ip - 1)*idot
         end do
         if (na == 0) return
         n2 = n + n
         do concurrent(i=1:n2)
            c(i) = ch(i)
         end do
      end subroutine cfftf1
