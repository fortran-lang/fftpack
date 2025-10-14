      subroutine cffti1(n, wa, ifac)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         integer, intent(out) :: ifac(*)
         real(dp), intent(out) :: wa(*)
         real(dp) :: arg, argh, argld, fi
         integer :: i, i1, ib, ido, idot, ii, ip, ipm, j, k1, &
                    l1, l2, ld, nf, nl, nq, nr, ntry
         integer, dimension(4), parameter :: ntryh = [3, 4, 2, 5]
         real(dp), parameter :: tpi = 2.0_dp*acos(-1.0_dp) ! 2 * pi
         nl = n
         nf = 0
         j = 0
100      j = j + 1
         if (j <= 4) then
            ntry = ntryh(j)
         else
            ntry = ntry + 2
         end if
200      nq = nl/ntry
         nr = nl - ntry*nq
         if (nr /= 0) goto 100
         nf = nf + 1
         ifac(nf + 2) = ntry
         nl = nq
         if (ntry == 2) then
            if (nf /= 1) then
               do i = 2, nf
                  ib = nf - i + 2
                  ifac(ib + 2) = ifac(ib + 1)
               end do
               ifac(3) = 2
            end if
         end if
         if (nl /= 1) goto 200
         ifac(1) = n
         ifac(2) = nf
         argh = tpi/real(n, kind=dp)
         i = 2
         l1 = 1
         do k1 = 1, nf
            ip = ifac(k1 + 2)
            ld = 0
            l2 = l1*ip
            ido = n/l2
            idot = ido + ido + 2
            ipm = ip - 1
            do j = 1, ipm
               i1 = i
               wa(i - 1) = 1.0_dp
               wa(i) = 0.0_dp
               ld = ld + l1
               fi = 0.0_dp
               argld = real(ld, kind=dp)*argh
               do ii = 4, idot, 2
                  i = i + 2
                  fi = fi + 1.0_dp
                  arg = fi*argld
                  wa(i - 1) = cos(arg)
                  wa(i) = sin(arg)
               end do
               if (ip > 5) then
                  wa(i1 - 1) = wa(i - 1)
                  wa(i1) = wa(i)
               end if
            end do
            l1 = l2
         end do
      end subroutine cffti1
