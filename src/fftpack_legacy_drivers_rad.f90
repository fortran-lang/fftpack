module fftpack_legacy_drivers_rad
   use fftpack_kinds, only: dp
   implicit none(type, external)
   private

   public :: radbg, radb2, radb3, radb4, radb5
   public :: radfg, radf2, radf3, radf4, radf5

   interface
      pure module subroutine radbg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
         implicit none(type, external)
         integer, intent(in) :: ido, ip, l1, idl1
         real(dp), intent(in) :: cc(ido, ip, l1), wa(*)
         real(dp), intent(inout) :: c1(ido, l1, ip), ch2(idl1, ip)
         real(dp), intent(out) :: c2(idl1, ip), ch(ido, l1, ip)
      end subroutine radbg
   end interface

   interface
      pure module subroutine radb2(ido, l1, cc, ch, wa1)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 2, l1), wa1(*)
         real(dp), intent(out) :: ch(ido, l1, 2)
      end subroutine radb2
   end interface

   interface
      pure module subroutine radb3(ido, l1, cc, ch, wa1, wa2)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 3, l1), wa1(*), wa2(*)
         real(dp), intent(out) :: ch(ido, l1, 3)
      end subroutine radb3
   end interface

   interface
      pure module subroutine radb4(ido, l1, cc, ch, wa1, wa2, wa3)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 4, l1), wa1(*), wa2(*), wa3(*)
         real(dp), intent(out) :: ch(ido, l1, 4)
      end subroutine radb4
   end interface

   interface
      pure module subroutine radb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 5, l1), wa1(*), wa2(*), wa3(*), wa4(*)
         real(dp), intent(out) :: ch(ido, l1, 5)
      end subroutine radb5
   end interface

   interface
      pure module subroutine radfg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
         implicit none(type, external)
         integer, intent(in) :: ido, ip, l1, idl1
         real(dp), intent(in) :: wa(*)
         real(dp), intent(inout) :: cc(ido, ip, l1)
         real(dp), intent(inout) :: c1(ido, l1, ip)
         real(dp), intent(inout) :: c2(idl1, ip)
         real(dp), intent(out) :: ch(ido, l1, ip)
         real(dp), intent(inout) :: ch2(idl1, ip)
      end subroutine radfg
   end interface

   interface
      pure module subroutine radf2(ido, l1, cc, ch, wa1)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, l1, 2), wa1(*)
         real(dp), intent(out) :: ch(ido, 2, l1)
      end subroutine radf2
   end interface

   interface
      pure module subroutine radf3(ido, l1, cc, ch, wa1, wa2)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, l1, 3), wa1(*), wa2(*)
         real(dp), intent(out) :: ch(ido, 3, l1)
      end subroutine radf3
   end interface

   interface
      pure module subroutine radf4(ido, l1, cc, ch, wa1, wa2, wa3)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, l1, 4), wa1(*), wa2(*), wa3(*)
         real(dp), intent(out) :: ch(ido, 4, l1)
      end subroutine radf4
   end interface

   interface
      pure module subroutine radf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, l1, 5), wa1(*), wa2(*), wa3(*), wa4(*)
         real(dp), intent(out) :: ch(ido, 5, l1)
      end subroutine radf5
   end interface

contains

   pure module subroutine radbg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
      implicit none(type, external)
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
         do concurrent(k=1:l1, i=1:ido)
            ch(i, k, 1) = cc(i, 1, k)
         end do
      else
         do concurrent(k=1:l1, i=1:ido)
            ch(i, k, 1) = cc(i, 1, k)
         end do
      end if
      do concurrent(k=1:l1, j=2:ipph)
         jc = ipp2 - j
         j2 = j + j
         ch(1, k, j) = cc(ido, j2 - 2, k) + cc(ido, j2 - 2, k)
         ch(1, k, jc) = cc(1, j2 - 1, k) + cc(1, j2 - 1, k)
      end do
      if (ido /= 1) then
         if (nbd < l1) then
            do concurrent(k=1:l1, j=2:ipph, i=3:ido:2)
               jc = ipp2 - j
               ic = idp2 - i
               ch(i - 1, k, j) = cc(i - 1, 2*j - 1, k) + cc(ic - 1, 2*j - 2, k)
               ch(i - 1, k, jc) = cc(i - 1, 2*j - 1, k) - cc(ic - 1, 2*j - 2, k)
               ch(i, k, j) = cc(i, 2*j - 1, k) - cc(ic, 2*j - 2, k)
               ch(i, k, jc) = cc(i, 2*j - 1, k) + cc(ic, 2*j - 2, k)
            end do
         else
            do concurrent(k=1:l1, j=2:ipph, i=3:ido:2)
               jc = ipp2 - j
               ic = idp2 - i
               ch(i - 1, k, j) = cc(i - 1, 2*j - 1, k) + cc(ic - 1, 2*j - 2, k)
               ch(i - 1, k, jc) = cc(i - 1, 2*j - 1, k) - cc(ic - 1, 2*j - 2, k)
               ch(i, k, j) = cc(i, 2*j - 1, k) - cc(ic, 2*j - 2, k)
               ch(i, k, jc) = cc(i, 2*j - 1, k) + cc(ic, 2*j - 2, k)
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
         do concurrent(ik=1:idl1)
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
            do concurrent(ik=1:idl1)
               c2(ik, l) = c2(ik, l) + ar2*ch2(ik, j)
               c2(ik, lc) = c2(ik, lc) + ai2*ch2(ik, jc)
            end do
         end do
      end do
      do concurrent(ik=1:idl1, j=2:ipph)
         ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
      end do
      do concurrent(j=2:ipph, k=1:l1)
         jc = ipp2 - j
         ch(1, k, j) = c1(1, k, j) - c1(1, k, jc)
         ch(1, k, jc) = c1(1, k, j) + c1(1, k, jc)
      end do
      if (ido /= 1) then
         if (nbd < l1) then
            do concurrent(j=2:ipph, k=1:l1, i=3:ido:2)
               jc = ipp2 - j
               ch(i - 1, k, j) = c1(i - 1, k, j) - c1(i, k, jc)
               ch(i - 1, k, jc) = c1(i - 1, k, j) + c1(i, k, jc)
               ch(i, k, j) = c1(i, k, j) + c1(i - 1, k, jc)
               ch(i, k, jc) = c1(i, k, j) - c1(i - 1, k, jc)
            end do
         else
            do concurrent(j=2:ipph, k=1:l1, i=3:ido:2)
               jc = ipp2 - j
               ch(i - 1, k, j) = c1(i - 1, k, j) - c1(i, k, jc)
               ch(i - 1, k, jc) = c1(i - 1, k, j) + c1(i, k, jc)
               ch(i, k, j) = c1(i, k, j) + c1(i - 1, k, jc)
               ch(i, k, jc) = c1(i, k, j) - c1(i - 1, k, jc)
            end do
         end if
      end if
      if (ido == 1) return
      do concurrent(ik=1:idl1)
         c2(ik, 1) = ch2(ik, 1)
      end do
      do concurrent(j=2:ip, k=1:l1)
         c1(1, k, j) = ch(1, k, j)
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

   pure module subroutine radb2(ido, l1, cc, ch, wa1)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 2, l1), wa1(*)
      real(dp), intent(out) :: ch(ido, l1, 2)
      real(dp) :: ti2, tr2
      integer :: i, ic, idp2, k
      do concurrent(k=1:l1)
         ch(1, k, 1) = cc(1, 1, k) + cc(ido, 2, k)
         ch(1, k, 2) = cc(1, 1, k) - cc(ido, 2, k)
      end do
      if (ido < 2) return
      if (ido /= 2) then
         idp2 = ido + 2
         do concurrent(k=1:l1, i=3:ido:2)
            ic = idp2 - i
            ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(ic - 1, 2, k)
            tr2 = cc(i - 1, 1, k) - cc(ic - 1, 2, k)
            ch(i, k, 1) = cc(i, 1, k) - cc(ic, 2, k)
            ti2 = cc(i, 1, k) + cc(ic, 2, k)
            ch(i - 1, k, 2) = wa1(i - 2)*tr2 - wa1(i - 1)*ti2
            ch(i, k, 2) = wa1(i - 2)*ti2 + wa1(i - 1)*tr2
         end do
         if (mod(ido, 2) == 1) return
      end if
      do concurrent(k=1:l1)
         ch(ido, k, 1) = cc(ido, 1, k) + cc(ido, 1, k)
         ch(ido, k, 2) = -(cc(1, 2, k) + cc(1, 2, k))
      end do
   end subroutine radb2

   pure module subroutine radb3(ido, l1, cc, ch, wa1, wa2)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 3, l1), wa1(*), wa2(*)
      real(dp), intent(out) :: ch(ido, l1, 3)
      real(dp) :: ci2, ci3, cr2, cr3, di2, di3, &
                  dr2, dr3, ti2, tr2
      integer :: i, ic, idp2, k
      real(dp), parameter :: taur = -0.5_dp
      real(dp), parameter :: taui = sqrt(3.0_dp)/2.0_dp
      do concurrent(k=1:l1)
         tr2 = cc(ido, 2, k) + cc(ido, 2, k)
         cr2 = cc(1, 1, k) + taur*tr2
         ch(1, k, 1) = cc(1, 1, k) + tr2
         ci3 = taui*(cc(1, 3, k) + cc(1, 3, k))
         ch(1, k, 2) = cr2 - ci3
         ch(1, k, 3) = cr2 + ci3
      end do
      if (ido == 1) return
      idp2 = ido + 2
      do concurrent(k=1:l1, i=3:ido:2)
         ic = idp2 - i
         tr2 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
         cr2 = cc(i - 1, 1, k) + taur*tr2
         ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
         ti2 = cc(i, 3, k) - cc(ic, 2, k)
         ci2 = cc(i, 1, k) + taur*ti2
         ch(i, k, 1) = cc(i, 1, k) + ti2
         cr3 = taui*(cc(i - 1, 3, k) - cc(ic - 1, 2, k))
         ci3 = taui*(cc(i, 3, k) + cc(ic, 2, k))
         dr2 = cr2 - ci3
         dr3 = cr2 + ci3
         di2 = ci2 + cr3
         di3 = ci2 - cr3
         ch(i - 1, k, 2) = wa1(i - 2)*dr2 - wa1(i - 1)*di2
         ch(i, k, 2) = wa1(i - 2)*di2 + wa1(i - 1)*dr2
         ch(i - 1, k, 3) = wa2(i - 2)*dr3 - wa2(i - 1)*di3
         ch(i, k, 3) = wa2(i - 2)*di3 + wa2(i - 1)*dr3
      end do
   end subroutine radb3

   pure module subroutine radb4(ido, l1, cc, ch, wa1, wa2, wa3)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 4, l1), wa1(*), wa2(*), wa3(*)
      real(dp), intent(out) :: ch(ido, l1, 4)
      real(dp) :: ci2, ci3, ci4, cr2, cr3, cr4, &
                  ti1, ti2, ti3, ti4, tr1, tr2, tr3, &
                  tr4
      integer :: i, ic, idp2, k
      real(dp), parameter :: sqrt2 = sqrt(2.0_dp)
      do concurrent(k=1:l1)
         tr1 = cc(1, 1, k) - cc(ido, 4, k)
         tr2 = cc(1, 1, k) + cc(ido, 4, k)
         tr3 = cc(ido, 2, k) + cc(ido, 2, k)
         tr4 = cc(1, 3, k) + cc(1, 3, k)
         ch(1, k, 1) = tr2 + tr3
         ch(1, k, 2) = tr1 - tr4
         ch(1, k, 3) = tr2 - tr3
         ch(1, k, 4) = tr1 + tr4
      end do
      if (ido < 2) return
      if (ido /= 2) then
         idp2 = ido + 2
         do concurrent(k=1:l1, i=3:ido:2)
            ic = idp2 - i
            ti1 = cc(i, 1, k) + cc(ic, 4, k)
            ti2 = cc(i, 1, k) - cc(ic, 4, k)
            ti3 = cc(i, 3, k) - cc(ic, 2, k)
            tr4 = cc(i, 3, k) + cc(ic, 2, k)
            tr1 = cc(i - 1, 1, k) - cc(ic - 1, 4, k)
            tr2 = cc(i - 1, 1, k) + cc(ic - 1, 4, k)
            ti4 = cc(i - 1, 3, k) - cc(ic - 1, 2, k)
            tr3 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
            ch(i - 1, k, 1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i, k, 1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 - tr4
            cr4 = tr1 + tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i - 1, k, 2) = wa1(i - 2)*cr2 - wa1(i - 1)*ci2
            ch(i, k, 2) = wa1(i - 2)*ci2 + wa1(i - 1)*cr2
            ch(i - 1, k, 3) = wa2(i - 2)*cr3 - wa2(i - 1)*ci3
            ch(i, k, 3) = wa2(i - 2)*ci3 + wa2(i - 1)*cr3
            ch(i - 1, k, 4) = wa3(i - 2)*cr4 - wa3(i - 1)*ci4
            ch(i, k, 4) = wa3(i - 2)*ci4 + wa3(i - 1)*cr4
         end do
         if (mod(ido, 2) == 1) return
      end if
      do concurrent(k=1:l1)
         ti1 = cc(1, 2, k) + cc(1, 4, k)
         ti2 = cc(1, 4, k) - cc(1, 2, k)
         tr1 = cc(ido, 1, k) - cc(ido, 3, k)
         tr2 = cc(ido, 1, k) + cc(ido, 3, k)
         ch(ido, k, 1) = tr2 + tr2
         ch(ido, k, 2) = sqrt2*(tr1 - ti1)
         ch(ido, k, 3) = ti2 + ti2
         ch(ido, k, 4) = -sqrt2*(tr1 + ti1)
      end do
   end subroutine radb4

   pure module subroutine radb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 5, l1), wa1(*), wa2(*), wa3(*), wa4(*)
      real(dp), intent(out) :: ch(ido, l1, 5)
      real(dp) :: ci2, ci3, ci4, ci5, cr2, cr3, &
                  cr4, cr5, di2, di3, di4, di5, dr2, dr3, &
                  dr4, dr5
      real(dp) :: ti2, ti3, ti4, ti5, tr2, tr3, &
                  tr4, tr5
      integer :: i, ic, idp2, k
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp), parameter :: tr11 = cos(2.0_dp*pi/5.0_dp)
      real(dp), parameter :: ti11 = sin(2.0_dp*pi/5.0_dp)
      real(dp), parameter :: tr12 = cos(4.0_dp*pi/5.0_dp)
      real(dp), parameter :: ti12 = sin(4.0_dp*pi/5.0_dp)
      do concurrent(k=1:l1)
         ti5 = cc(1, 3, k) + cc(1, 3, k)
         ti4 = cc(1, 5, k) + cc(1, 5, k)
         tr2 = cc(ido, 2, k) + cc(ido, 2, k)
         tr3 = cc(ido, 4, k) + cc(ido, 4, k)
         ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
         cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
         cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
         ci5 = ti11*ti5 + ti12*ti4
         ci4 = ti12*ti5 - ti11*ti4
         ch(1, k, 2) = cr2 - ci5
         ch(1, k, 3) = cr3 - ci4
         ch(1, k, 4) = cr3 + ci4
         ch(1, k, 5) = cr2 + ci5
      end do
      if (ido == 1) return
      idp2 = ido + 2
      do concurrent(k=1:l1, i=3:ido:2)
         ic = idp2 - i
         ti5 = cc(i, 3, k) + cc(ic, 2, k)
         ti2 = cc(i, 3, k) - cc(ic, 2, k)
         ti4 = cc(i, 5, k) + cc(ic, 4, k)
         ti3 = cc(i, 5, k) - cc(ic, 4, k)
         tr5 = cc(i - 1, 3, k) - cc(ic - 1, 2, k)
         tr2 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
         tr4 = cc(i - 1, 5, k) - cc(ic - 1, 4, k)
         tr3 = cc(i - 1, 5, k) + cc(ic - 1, 4, k)
         ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2 + tr3
         ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
         cr2 = cc(i - 1, 1, k) + tr11*tr2 + tr12*tr3
         ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
         cr3 = cc(i - 1, 1, k) + tr12*tr2 + tr11*tr3
         ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
         cr5 = ti11*tr5 + ti12*tr4
         ci5 = ti11*ti5 + ti12*ti4
         cr4 = ti12*tr5 - ti11*tr4
         ci4 = ti12*ti5 - ti11*ti4
         dr3 = cr3 - ci4
         dr4 = cr3 + ci4
         di3 = ci3 + cr4
         di4 = ci3 - cr4
         dr5 = cr2 + ci5
         dr2 = cr2 - ci5
         di5 = ci2 - cr5
         di2 = ci2 + cr5
         ch(i - 1, k, 2) = wa1(i - 2)*dr2 - wa1(i - 1)*di2
         ch(i, k, 2) = wa1(i - 2)*di2 + wa1(i - 1)*dr2
         ch(i - 1, k, 3) = wa2(i - 2)*dr3 - wa2(i - 1)*di3
         ch(i, k, 3) = wa2(i - 2)*di3 + wa2(i - 1)*dr3
         ch(i - 1, k, 4) = wa3(i - 2)*dr4 - wa3(i - 1)*di4
         ch(i, k, 4) = wa3(i - 2)*di4 + wa3(i - 1)*dr4
         ch(i - 1, k, 5) = wa4(i - 2)*dr5 - wa4(i - 1)*di5
         ch(i, k, 5) = wa4(i - 2)*di5 + wa4(i - 1)*dr5
      end do
   end subroutine radb5

   pure module subroutine radfg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
      implicit none(type, external)
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

   pure module subroutine radf2(ido, l1, cc, ch, wa1)
      implicit none(type, external)
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
         do concurrent(k=1:l1, i=3:ido:2)
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

   pure module subroutine radf3(ido, l1, cc, ch, wa1, wa2)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, l1, 3), wa1(*), wa2(*)
      real(dp), intent(out) :: ch(ido, 3, l1)
      real(dp) :: ci2, cr2, di2, di3, dr2, dr3, &
                  ti2, ti3, tr2, tr3
      integer :: i, ic, idp2, k
      real(dp), parameter :: taur = -0.5_dp
      ! note: original comment said this was -sqrt(3)/2 but value was 0.86602540378443864676d0
      real(dp), parameter :: taui = sqrt(3.0_dp)/2.0_dp
      do concurrent(k=1:l1)
         cr2 = cc(1, k, 2) + cc(1, k, 3)
         ch(1, 1, k) = cc(1, k, 1) + cr2
         ch(1, 3, k) = taui*(cc(1, k, 3) - cc(1, k, 2))
         ch(ido, 2, k) = cc(1, k, 1) + taur*cr2
      end do
      if (ido == 1) return
      idp2 = ido + 2
      do concurrent(k=1:l1, i=3:ido:2)
         ic = idp2 - i
         dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
         di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
         dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
         di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
         cr2 = dr2 + dr3
         ci2 = di2 + di3
         ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2
         ch(i, 1, k) = cc(i, k, 1) + ci2
         tr2 = cc(i - 1, k, 1) + taur*cr2
         ti2 = cc(i, k, 1) + taur*ci2
         tr3 = taui*(di2 - di3)
         ti3 = taui*(dr3 - dr2)
         ch(i - 1, 3, k) = tr2 + tr3
         ch(ic - 1, 2, k) = tr2 - tr3
         ch(i, 3, k) = ti2 + ti3
         ch(ic, 2, k) = ti3 - ti2
      end do
   end subroutine radf3

   pure module subroutine radf4(ido, l1, cc, ch, wa1, wa2, wa3)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, l1, 4), wa1(*), wa2(*), wa3(*)
      real(dp), intent(out) :: ch(ido, 4, l1)
      real(dp) :: ci2, ci3, ci4, cr2, cr3, cr4, &
                  ti1, ti2, ti3, ti4, tr1, tr2, tr3, &
                  tr4
      integer :: i, ic, idp2, k
      real(dp), parameter :: hsqt2 = sqrt(2.0_dp)/2.0_dp
      do concurrent(k=1:l1)
         tr1 = cc(1, k, 2) + cc(1, k, 4)
         tr2 = cc(1, k, 1) + cc(1, k, 3)
         ch(1, 1, k) = tr1 + tr2
         ch(ido, 4, k) = tr2 - tr1
         ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 3)
         ch(1, 3, k) = cc(1, k, 4) - cc(1, k, 2)
      end do
      if (ido < 2) return
      if (ido /= 2) then
         idp2 = ido + 2
         do concurrent(k=1:l1, i=3:ido:2)
            ic = idp2 - i
            cr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
            ci2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
            cr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
            ci3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
            cr4 = wa3(i - 2)*cc(i - 1, k, 4) + wa3(i - 1)*cc(i, k, 4)
            ci4 = wa3(i - 2)*cc(i, k, 4) - wa3(i - 1)*cc(i - 1, k, 4)
            tr1 = cr2 + cr4
            tr4 = cr4 - cr2
            ti1 = ci2 + ci4
            ti4 = ci2 - ci4
            ti2 = cc(i, k, 1) + ci3
            ti3 = cc(i, k, 1) - ci3
            tr2 = cc(i - 1, k, 1) + cr3
            tr3 = cc(i - 1, k, 1) - cr3
            ch(i - 1, 1, k) = tr1 + tr2
            ch(ic - 1, 4, k) = tr2 - tr1
            ch(i, 1, k) = ti1 + ti2
            ch(ic, 4, k) = ti1 - ti2
            ch(i - 1, 3, k) = ti4 + tr3
            ch(ic - 1, 2, k) = tr3 - ti4
            ch(i, 3, k) = tr4 + ti3
            ch(ic, 2, k) = tr4 - ti3
         end do
         if (mod(ido, 2) == 1) return
      end if
      do concurrent(k=1:l1)
         ti1 = -hsqt2*(cc(ido, k, 2) + cc(ido, k, 4))
         tr1 = hsqt2*(cc(ido, k, 2) - cc(ido, k, 4))
         ch(ido, 1, k) = tr1 + cc(ido, k, 1)
         ch(ido, 3, k) = cc(ido, k, 1) - tr1
         ch(1, 2, k) = ti1 - cc(ido, k, 3)
         ch(1, 4, k) = ti1 + cc(ido, k, 3)
      end do
   end subroutine radf4

   pure module subroutine radf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, l1, 5), wa1(*), wa2(*), wa3(*), wa4(*)
      real(dp), intent(out) :: ch(ido, 5, l1)
      real(dp) :: ci2, ci3, ci4, ci5, cr2, cr3, &
                  cr4, cr5, di2, di3, di4, di5, dr2, dr3, &
                  dr4, dr5
      real(dp) :: ti2, ti3, ti4, ti5, tr2, tr3, &
                  tr4, tr5
      integer :: i, ic, idp2, k
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp), parameter :: tr11 = cos(2.0_dp*pi/5.0_dp)
      real(dp), parameter :: ti11 = sin(2.0_dp*pi/5.0_dp)
      real(dp), parameter :: tr12 = cos(4.0_dp*pi/5.0_dp)
      real(dp), parameter :: ti12 = sin(4.0_dp*pi/5.0_dp)
      do concurrent(k=1:l1)
         cr2 = cc(1, k, 5) + cc(1, k, 2)
         ci5 = cc(1, k, 5) - cc(1, k, 2)
         cr3 = cc(1, k, 4) + cc(1, k, 3)
         ci4 = cc(1, k, 4) - cc(1, k, 3)
         ch(1, 1, k) = cc(1, k, 1) + cr2 + cr3
         ch(ido, 2, k) = cc(1, k, 1) + tr11*cr2 + tr12*cr3
         ch(1, 3, k) = ti11*ci5 + ti12*ci4
         ch(ido, 4, k) = cc(1, k, 1) + tr12*cr2 + tr11*cr3
         ch(1, 5, k) = ti12*ci5 - ti11*ci4
      end do
      if (ido == 1) return
      idp2 = ido + 2
      do concurrent(k=1:l1, i=3:ido:2)
         ic = idp2 - i
         dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
         di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
         dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
         di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
         dr4 = wa3(i - 2)*cc(i - 1, k, 4) + wa3(i - 1)*cc(i, k, 4)
         di4 = wa3(i - 2)*cc(i, k, 4) - wa3(i - 1)*cc(i - 1, k, 4)
         dr5 = wa4(i - 2)*cc(i - 1, k, 5) + wa4(i - 1)*cc(i, k, 5)
         di5 = wa4(i - 2)*cc(i, k, 5) - wa4(i - 1)*cc(i - 1, k, 5)
         cr2 = dr2 + dr5
         ci5 = dr5 - dr2
         cr5 = di2 - di5
         ci2 = di2 + di5
         cr3 = dr3 + dr4
         ci4 = dr4 - dr3
         cr4 = di3 - di4
         ci3 = di3 + di4
         ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2 + cr3
         ch(i, 1, k) = cc(i, k, 1) + ci2 + ci3
         tr2 = cc(i - 1, k, 1) + tr11*cr2 + tr12*cr3
         ti2 = cc(i, k, 1) + tr11*ci2 + tr12*ci3
         tr3 = cc(i - 1, k, 1) + tr12*cr2 + tr11*cr3
         ti3 = cc(i, k, 1) + tr12*ci2 + tr11*ci3
         tr5 = ti11*cr5 + ti12*cr4
         ti5 = ti11*ci5 + ti12*ci4
         tr4 = ti12*cr5 - ti11*cr4
         ti4 = ti12*ci5 - ti11*ci4
         ch(i - 1, 3, k) = tr2 + tr5
         ch(ic - 1, 2, k) = tr2 - tr5
         ch(i, 3, k) = ti2 + ti5
         ch(ic, 2, k) = ti5 - ti2
         ch(i - 1, 5, k) = tr3 + tr4
         ch(ic - 1, 4, k) = tr3 - tr4
         ch(i, 5, k) = ti3 + ti4
         ch(ic, 4, k) = ti4 - ti3
      end do
   end subroutine radf5
end module fftpack_legacy_drivers_rad
