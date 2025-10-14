      subroutine rffti1(n, wa, ifac)
         use fftpack_kind, only: dp => rk
         implicit none
         integer, intent(in) :: n
         real(dp), intent(out) :: wa(*)
         integer, intent(out) :: ifac(*)
         real(dp) :: arg, argh, argld, fi
         integer :: i, ib, ido, ii, ip, ipm, is, j, k1, l1, &
                    l2, ld, nf, nfm1, nl, nq, nr, ntry
         integer, dimension(4), parameter :: ntryh = [4, 2, 3, 5]
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
         argh = tpi/real(n, dp)
         is = 0
         nfm1 = nf - 1
         l1 = 1
         if (nfm1 == 0) return
         do k1 = 1, nfm1
            ip = ifac(k1 + 2)
            ld = 0
            l2 = l1*ip
            ido = n/l2
            ipm = ip - 1
            do j = 1, ipm
               ld = ld + l1
               i = is
               argld = real(ld, dp)*argh
               fi = 0.0_dp
               do ii = 3, ido, 2
                  i = i + 2
                  fi = fi + 1.0_dp
                  arg = fi*argld
                  wa(i - 1) = cos(arg)
                  wa(i) = sin(arg)
               end do
               is = is + ido
            end do
            l1 = l2
         end do
      end subroutine rffti1
