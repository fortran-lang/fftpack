      subroutine rfftb1(n, c, ch, wa, ifac)
         use fftpack_kinds, only: dp
         use fftpack_legacy_drivers_rad, only: radbg, radb2, radb3, radb4, radb5
         implicit none
         integer, intent(in) :: n
         real(dp), intent(inout) :: c(*)
         real(dp), intent(in) :: wa(*)
         real(dp), intent(inout) :: ch(*)
         integer, intent(in) :: ifac(*)
         integer :: i, idl1, ido, ip, iw, ix2, ix3, ix4, k1, &
                    l1, l2, na, nf
         nf = ifac(2)
         na = 0
         l1 = 1
         iw = 1
         do k1 = 1, nf
            ip = ifac(k1 + 2)
            l2 = ip*l1
            ido = n/l2
            idl1 = ido*l1
            if (ip == 4) then
               ix2 = iw + ido
               ix3 = ix2 + ido
               if (na /= 0) then
                  call radb4(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
               else
                  call radb4(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
               end if
               na = 1 - na
            elseif (ip == 2) then
               if (na /= 0) then
                  call radb2(ido, l1, ch, c, wa(iw))
               else
                  call radb2(ido, l1, c, ch, wa(iw))
               end if
               na = 1 - na
            elseif (ip == 3) then
               ix2 = iw + ido
               if (na /= 0) then
                  call radb3(ido, l1, ch, c, wa(iw), wa(ix2))
               else
                  call radb3(ido, l1, c, ch, wa(iw), wa(ix2))
               end if
               na = 1 - na
            elseif (ip /= 5) then
               if (na /= 0) then
                  call radbg(ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
               else
                  call radbg(ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
               end if
               if (ido == 1) na = 1 - na
            else
               ix2 = iw + ido
               ix3 = ix2 + ido
               ix4 = ix3 + ido
               if (na /= 0) then
                  call radb5(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
               else
                  call radb5(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
               end if
               na = 1 - na
            end if
            l1 = l2
            iw = iw + (ip - 1)*ido
         end do
         if (na == 0) return
         do concurrent(i=1:n)
            c(i) = ch(i)
         end do
      end subroutine rfftb1
