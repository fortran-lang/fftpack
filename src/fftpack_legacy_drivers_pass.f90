module fftpack_legacy_drivers_pass
   use fftpack_kinds, only: dp
   implicit none(type, external)
   private

   public :: passb, passb2, passb3, passb4, passb5
   public :: passf, passf2, passf3, passf4, passf5

   interface
      pure module subroutine passb(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
         implicit none(type, external)
         integer, intent(out) :: nac
         integer, intent(in) :: ido, ip, l1, idl1
         real(dp), intent(in) :: cc(ido, ip, l1), wa(*)
         real(dp), intent(out) :: c1(ido, l1, ip), c2(idl1, ip), ch(ido, l1, ip)
         real(dp), intent(inout) :: ch2(idl1, ip)
      end subroutine passb
   end interface

   interface
      pure module subroutine passb2(ido, l1, cc, ch, wa1)
         implicit none
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 2, l1), wa1(*)
         real(dp), intent(out) :: ch(ido, l1, 2)
      end subroutine passb2
   end interface

   interface
      pure module subroutine passb3(ido, l1, cc, ch, wa1, wa2)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 3, l1), wa1(*), wa2(*)
         real(dp), intent(out) :: ch(ido, l1, 3)
      end subroutine passb3
   end interface

   interface
      pure module subroutine passb4(ido, l1, cc, ch, wa1, wa2, wa3)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 4, l1), wa1(*), wa2(*), wa3(*)
         real(dp), intent(out) :: ch(ido, l1, 4)
      end subroutine passb4
   end interface

   interface
      pure module subroutine passb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 5, l1), wa1(*), wa2(*), wa3(*), wa4(*)
         real(dp), intent(out) :: ch(ido, l1, 5)
      end subroutine passb5
   end interface

   interface
      pure module subroutine passf(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
         implicit none(type, external)
         integer, intent(out) :: nac
         integer, intent(in) :: ido, ip, l1, idl1
         real(dp), intent(in) :: cc(ido, ip, l1), wa(*)
         real(dp), intent(out) :: c1(ido, l1, ip), c2(idl1, ip), ch(ido, l1, ip)
         real(dp), intent(inout) :: ch2(idl1, ip)
      end subroutine passf
   end interface

   interface
      pure module subroutine passf2(ido, l1, cc, ch, wa1)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 2, l1), wa1(*)
         real(dp), intent(out) :: ch(ido, l1, 2)
      end subroutine passf2
   end interface

   interface
      pure module subroutine passf3(ido, l1, cc, ch, wa1, wa2)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 3, l1), wa1(*), wa2(*)
         real(dp), intent(out) :: ch(ido, l1, 3)
      end subroutine passf3
   end interface

   interface
      pure module subroutine passf4(ido, l1, cc, ch, wa1, wa2, wa3)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 4, l1), wa1(*), wa2(*), wa3(*)
         real(dp), intent(out) :: ch(ido, l1, 4)
      end subroutine passf4
   end interface

   interface
      pure module subroutine passf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
         implicit none(type, external)
         integer, intent(in) :: ido, l1
         real(dp), intent(in) :: cc(ido, 5, l1), wa1(*), wa2(*), wa3(*), wa4(*)
         real(dp), intent(out) :: ch(ido, l1, 5)
      end subroutine passf5
   end interface

contains

   pure module subroutine passb(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
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
            c2(ik, lc) = wa(idl)*ch2(ik, ip)
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
               c2(ik, lc) = c2(ik, lc) + wai*ch2(ik, jc)
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
                  c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij) &
                                    *ch(i, k, j)
                  c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij) &
                                *ch(i - 1, k, j)
               end do
            end do
         end do
         return
      end if
      idij = 0
      do j = 2, ip
         idij = idij + 2
         do i = 4, ido, 2
            idij = idij + 2
            do concurrent(k=1:l1)
               c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij)*ch(i, k, j)
               c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
            end do
         end do
      end do
   end subroutine passb

   pure module subroutine passb2(ido, l1, cc, ch, wa1)
      implicit none(type, external)
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
            ch(i, k, 2) = wa1(i - 1)*ti2 + wa1(i)*tr2
            ch(i - 1, k, 2) = wa1(i - 1)*tr2 - wa1(i)*ti2
         end do
      else
         do concurrent(k=1:l1)
            ch(1, k, 1) = cc(1, 1, k) + cc(1, 2, k)
            ch(1, k, 2) = cc(1, 1, k) - cc(1, 2, k)
            ch(2, k, 1) = cc(2, 1, k) + cc(2, 2, k)
            ch(2, k, 2) = cc(2, 1, k) - cc(2, 2, k)
         end do
      end if
   end subroutine passb2

   pure module subroutine passb3(ido, l1, cc, ch, wa1, wa2)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 3, l1), wa1(*), wa2(*)
      real(dp), intent(out) :: ch(ido, l1, 3)
      real(dp) :: ci2, ci3, cr2, cr3, di2, di3, &
                  dr2, dr3, ti2, tr2
      integer :: i, k
      real(dp), parameter :: taur = -0.5_dp
      real(dp), parameter :: taui = sqrt(3.0_dp)/2.0_dp
      if (ido /= 2) then
         do concurrent(k=1:l1, i=2:ido:2)
            tr2 = cc(i - 1, 2, k) + cc(i - 1, 3, k)
            cr2 = cc(i - 1, 1, k) + taur*tr2
            ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
            ti2 = cc(i, 2, k) + cc(i, 3, k)
            ci2 = cc(i, 1, k) + taur*ti2
            ch(i, k, 1) = cc(i, 1, k) + ti2
            cr3 = taui*(cc(i - 1, 2, k) - cc(i - 1, 3, k))
            ci3 = taui*(cc(i, 2, k) - cc(i, 3, k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i, k, 2) = wa1(i - 1)*di2 + wa1(i)*dr2
            ch(i - 1, k, 2) = wa1(i - 1)*dr2 - wa1(i)*di2
            ch(i, k, 3) = wa2(i - 1)*di3 + wa2(i)*dr3
            ch(i - 1, k, 3) = wa2(i - 1)*dr3 - wa2(i)*di3
         end do
      else
         do concurrent(k=1:l1)
            tr2 = cc(1, 2, k) + cc(1, 3, k)
            cr2 = cc(1, 1, k) + taur*tr2
            ch(1, k, 1) = cc(1, 1, k) + tr2
            ti2 = cc(2, 2, k) + cc(2, 3, k)
            ci2 = cc(2, 1, k) + taur*ti2
            ch(2, k, 1) = cc(2, 1, k) + ti2
            cr3 = taui*(cc(1, 2, k) - cc(1, 3, k))
            ci3 = taui*(cc(2, 2, k) - cc(2, 3, k))
            ch(1, k, 2) = cr2 - ci3
            ch(1, k, 3) = cr2 + ci3
            ch(2, k, 2) = ci2 + cr3
            ch(2, k, 3) = ci2 - cr3
         end do
      end if
   end subroutine passb3

   pure module subroutine passb4(ido, l1, cc, ch, wa1, wa2, wa3)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 4, l1), wa1(*), wa2(*), wa3(*)
      real(dp), intent(out) :: ch(ido, l1, 4)
      real(dp) :: ci2, ci3, ci4, cr2, cr3, cr4, &
                     & ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
      integer :: i, k
      if (ido /= 2) then
         do concurrent(k=1:l1, i=2:ido:2)
            ti1 = cc(i, 1, k) - cc(i, 3, k)
            ti2 = cc(i, 1, k) + cc(i, 3, k)
            ti3 = cc(i, 2, k) + cc(i, 4, k)
            tr4 = cc(i, 4, k) - cc(i, 2, k)
            tr1 = cc(i - 1, 1, k) - cc(i - 1, 3, k)
            tr2 = cc(i - 1, 1, k) + cc(i - 1, 3, k)
            ti4 = cc(i - 1, 2, k) - cc(i - 1, 4, k)
            tr3 = cc(i - 1, 2, k) + cc(i - 1, 4, k)
            ch(i - 1, k, 1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i, k, 1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 + tr4
            cr4 = tr1 - tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i - 1, k, 2) = wa1(i - 1)*cr2 - wa1(i)*ci2
            ch(i, k, 2) = wa1(i - 1)*ci2 + wa1(i)*cr2
            ch(i - 1, k, 3) = wa2(i - 1)*cr3 - wa2(i)*ci3
            ch(i, k, 3) = wa2(i - 1)*ci3 + wa2(i)*cr3
            ch(i - 1, k, 4) = wa3(i - 1)*cr4 - wa3(i)*ci4
            ch(i, k, 4) = wa3(i - 1)*ci4 + wa3(i)*cr4
         end do
      else
         do concurrent(k=1:l1)
            ti1 = cc(2, 1, k) - cc(2, 3, k)
            ti2 = cc(2, 1, k) + cc(2, 3, k)
            tr4 = cc(2, 4, k) - cc(2, 2, k)
            ti3 = cc(2, 2, k) + cc(2, 4, k)
            tr1 = cc(1, 1, k) - cc(1, 3, k)
            tr2 = cc(1, 1, k) + cc(1, 3, k)
            ti4 = cc(1, 2, k) - cc(1, 4, k)
            tr3 = cc(1, 2, k) + cc(1, 4, k)
            ch(1, k, 1) = tr2 + tr3
            ch(1, k, 3) = tr2 - tr3
            ch(2, k, 1) = ti2 + ti3
            ch(2, k, 3) = ti2 - ti3
            ch(1, k, 2) = tr1 + tr4
            ch(1, k, 4) = tr1 - tr4
            ch(2, k, 2) = ti1 + ti4
            ch(2, k, 4) = ti1 - ti4
         end do
      end if
   end subroutine passb4

   pure module subroutine passb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 5, l1), wa1(*), wa2(*), wa3(*), wa4(*)
      real(dp), intent(out) :: ch(ido, l1, 5)
      real(dp) :: ci2, ci3, ci4, ci5, cr2, cr3, &
                  cr4, cr5, di2, di3, di4, di5, dr2, dr3, &
                  dr4, dr5
      real(dp) :: ti2, ti3, ti4, ti5, tr2, tr3, &
                  tr4, tr5
      integer :: i, k
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp), parameter :: tr11 = cos(2.0_dp*pi/5.0_dp)
      real(dp), parameter :: ti11 = sin(2.0_dp*pi/5.0_dp)
      real(dp), parameter :: tr12 = cos(4.0_dp*pi/5.0_dp)
      real(dp), parameter :: ti12 = sin(4.0_dp*pi/5.0_dp)
      if (ido /= 2) then
         do concurrent(k=1:l1, i=2:ido:2)
            ti5 = cc(i, 2, k) - cc(i, 5, k)
            ti2 = cc(i, 2, k) + cc(i, 5, k)
            ti4 = cc(i, 3, k) - cc(i, 4, k)
            ti3 = cc(i, 3, k) + cc(i, 4, k)
            tr5 = cc(i - 1, 2, k) - cc(i - 1, 5, k)
            tr2 = cc(i - 1, 2, k) + cc(i - 1, 5, k)
            tr4 = cc(i - 1, 3, k) - cc(i - 1, 4, k)
            tr3 = cc(i - 1, 3, k) + cc(i - 1, 4, k)
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
            ch(i - 1, k, 2) = wa1(i - 1)*dr2 - wa1(i)*di2
            ch(i, k, 2) = wa1(i - 1)*di2 + wa1(i)*dr2
            ch(i - 1, k, 3) = wa2(i - 1)*dr3 - wa2(i)*di3
            ch(i, k, 3) = wa2(i - 1)*di3 + wa2(i)*dr3
            ch(i - 1, k, 4) = wa3(i - 1)*dr4 - wa3(i)*di4
            ch(i, k, 4) = wa3(i - 1)*di4 + wa3(i)*dr4
            ch(i - 1, k, 5) = wa4(i - 1)*dr5 - wa4(i)*di5
            ch(i, k, 5) = wa4(i - 1)*di5 + wa4(i)*dr5
         end do
      else
         do concurrent(k=1:l1)
            ti5 = cc(2, 2, k) - cc(2, 5, k)
            ti2 = cc(2, 2, k) + cc(2, 5, k)
            ti4 = cc(2, 3, k) - cc(2, 4, k)
            ti3 = cc(2, 3, k) + cc(2, 4, k)
            tr5 = cc(1, 2, k) - cc(1, 5, k)
            tr2 = cc(1, 2, k) + cc(1, 5, k)
            tr4 = cc(1, 3, k) - cc(1, 4, k)
            tr3 = cc(1, 3, k) + cc(1, 4, k)
            ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
            ch(2, k, 1) = cc(2, 1, k) + ti2 + ti3
            cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
            ci2 = cc(2, 1, k) + tr11*ti2 + tr12*ti3
            cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
            ci3 = cc(2, 1, k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            ch(1, k, 2) = cr2 - ci5
            ch(1, k, 5) = cr2 + ci5
            ch(2, k, 2) = ci2 + cr5
            ch(2, k, 3) = ci3 + cr4
            ch(1, k, 3) = cr3 - ci4
            ch(1, k, 4) = cr3 + ci4
            ch(2, k, 4) = ci3 - cr4
            ch(2, k, 5) = ci2 - cr5
         end do
      end if
   end subroutine passb5

   pure module subroutine passf(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
      implicit none(type, external)
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

   pure module subroutine passf2(ido, l1, cc, ch, wa1)
      implicit none(type, external)
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

   pure module subroutine passf3(ido, l1, cc, ch, wa1, wa2)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 3, l1), wa1(*), wa2(*)
      real(dp), intent(out) :: ch(ido, l1, 3)
      real(dp) :: ci2, ci3, cr2, cr3, di2, di3, &
                     & dr2, dr3, ti2, tr2
      integer :: i, k
      real(dp), parameter :: taur = -0.5_dp
      real(dp), parameter :: taui = -sqrt(3.0_dp)/2.0_dp
      if (ido /= 2) then
         do concurrent(k=1:l1, i=2:ido:2)
            tr2 = cc(i - 1, 2, k) + cc(i - 1, 3, k)
            cr2 = cc(i - 1, 1, k) + taur*tr2
            ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
            ti2 = cc(i, 2, k) + cc(i, 3, k)
            ci2 = cc(i, 1, k) + taur*ti2
            ch(i, k, 1) = cc(i, 1, k) + ti2
            cr3 = taui*(cc(i - 1, 2, k) - cc(i - 1, 3, k))
            ci3 = taui*(cc(i, 2, k) - cc(i, 3, k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i, k, 2) = wa1(i - 1)*di2 - wa1(i)*dr2
            ch(i - 1, k, 2) = wa1(i - 1)*dr2 + wa1(i)*di2
            ch(i, k, 3) = wa2(i - 1)*di3 - wa2(i)*dr3
            ch(i - 1, k, 3) = wa2(i - 1)*dr3 + wa2(i)*di3
         end do
      else
         do concurrent(k=1:l1)
            tr2 = cc(1, 2, k) + cc(1, 3, k)
            cr2 = cc(1, 1, k) + taur*tr2
            ch(1, k, 1) = cc(1, 1, k) + tr2
            ti2 = cc(2, 2, k) + cc(2, 3, k)
            ci2 = cc(2, 1, k) + taur*ti2
            ch(2, k, 1) = cc(2, 1, k) + ti2
            cr3 = taui*(cc(1, 2, k) - cc(1, 3, k))
            ci3 = taui*(cc(2, 2, k) - cc(2, 3, k))
            ch(1, k, 2) = cr2 - ci3
            ch(1, k, 3) = cr2 + ci3
            ch(2, k, 2) = ci2 + cr3
            ch(2, k, 3) = ci2 - cr3
         end do
      end if
   end subroutine passf3

   pure module subroutine passf4(ido, l1, cc, ch, wa1, wa2, wa3)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 4, l1), wa1(*), wa2(*), wa3(*)
      real(dp), intent(out) :: ch(ido, l1, 4)
      real(dp) :: ci2, ci3, ci4, cr2, cr3, cr4, &
                     & ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
      integer :: i, k
      if (ido /= 2) then
         do concurrent(k=1:l1, i=2:ido:2)
            ti1 = cc(i, 1, k) - cc(i, 3, k)
            ti2 = cc(i, 1, k) + cc(i, 3, k)
            ti3 = cc(i, 2, k) + cc(i, 4, k)
            tr4 = cc(i, 2, k) - cc(i, 4, k)
            tr1 = cc(i - 1, 1, k) - cc(i - 1, 3, k)
            tr2 = cc(i - 1, 1, k) + cc(i - 1, 3, k)
            ti4 = cc(i - 1, 4, k) - cc(i - 1, 2, k)
            tr3 = cc(i - 1, 2, k) + cc(i - 1, 4, k)
            ch(i - 1, k, 1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i, k, 1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 + tr4
            cr4 = tr1 - tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i - 1, k, 2) = wa1(i - 1)*cr2 + wa1(i)*ci2
            ch(i, k, 2) = wa1(i - 1)*ci2 - wa1(i)*cr2
            ch(i - 1, k, 3) = wa2(i - 1)*cr3 + wa2(i)*ci3
            ch(i, k, 3) = wa2(i - 1)*ci3 - wa2(i)*cr3
            ch(i - 1, k, 4) = wa3(i - 1)*cr4 + wa3(i)*ci4
            ch(i, k, 4) = wa3(i - 1)*ci4 - wa3(i)*cr4
         end do
      else
         do concurrent(k=1:l1)
            ti1 = cc(2, 1, k) - cc(2, 3, k)
            ti2 = cc(2, 1, k) + cc(2, 3, k)
            tr4 = cc(2, 2, k) - cc(2, 4, k)
            ti3 = cc(2, 2, k) + cc(2, 4, k)
            tr1 = cc(1, 1, k) - cc(1, 3, k)
            tr2 = cc(1, 1, k) + cc(1, 3, k)
            ti4 = cc(1, 4, k) - cc(1, 2, k)
            tr3 = cc(1, 2, k) + cc(1, 4, k)
            ch(1, k, 1) = tr2 + tr3
            ch(1, k, 3) = tr2 - tr3
            ch(2, k, 1) = ti2 + ti3
            ch(2, k, 3) = ti2 - ti3
            ch(1, k, 2) = tr1 + tr4
            ch(1, k, 4) = tr1 - tr4
            ch(2, k, 2) = ti1 + ti4
            ch(2, k, 4) = ti1 - ti4
         end do
      end if
   end subroutine passf4

   pure module subroutine passf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
      implicit none(type, external)
      integer, intent(in) :: ido, l1
      real(dp), intent(in) :: cc(ido, 5, l1), wa1(*), wa2(*), wa3(*), wa4(*)
      real(dp), intent(out) :: ch(ido, l1, 5)
      real(dp) :: ci2, ci3, ci4, ci5, cr2, cr3, &
                  cr4, cr5, di2, di3, di4, di5, dr2, dr3, &
                  dr4, dr5
      real(dp) :: ti2, ti3, ti4, ti5, tr2, tr3, &
                  tr4, tr5
      integer :: i, k
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp), parameter :: tr11 = cos(2.0_dp*pi/5.0_dp)
      real(dp), parameter :: ti11 = -sin(2.0_dp*pi/5.0_dp)
      real(dp), parameter :: tr12 = cos(4.0_dp*pi/5.0_dp)
      real(dp), parameter :: ti12 = -sin(4.0_dp*pi/5.0_dp)
      if (ido /= 2) then
         do concurrent(k=1:l1, i=2:ido:2)
            ti5 = cc(i, 2, k) - cc(i, 5, k)
            ti2 = cc(i, 2, k) + cc(i, 5, k)
            ti4 = cc(i, 3, k) - cc(i, 4, k)
            ti3 = cc(i, 3, k) + cc(i, 4, k)
            tr5 = cc(i - 1, 2, k) - cc(i - 1, 5, k)
            tr2 = cc(i - 1, 2, k) + cc(i - 1, 5, k)
            tr4 = cc(i - 1, 3, k) - cc(i - 1, 4, k)
            tr3 = cc(i - 1, 3, k) + cc(i - 1, 4, k)
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
            ch(i - 1, k, 2) = wa1(i - 1)*dr2 + wa1(i)*di2
            ch(i, k, 2) = wa1(i - 1)*di2 - wa1(i)*dr2
            ch(i - 1, k, 3) = wa2(i - 1)*dr3 + wa2(i)*di3
            ch(i, k, 3) = wa2(i - 1)*di3 - wa2(i)*dr3
            ch(i - 1, k, 4) = wa3(i - 1)*dr4 + wa3(i)*di4
            ch(i, k, 4) = wa3(i - 1)*di4 - wa3(i)*dr4
            ch(i - 1, k, 5) = wa4(i - 1)*dr5 + wa4(i)*di5
            ch(i, k, 5) = wa4(i - 1)*di5 - wa4(i)*dr5
         end do
      else
         do concurrent(k=1:l1)
            ti5 = cc(2, 2, k) - cc(2, 5, k)
            ti2 = cc(2, 2, k) + cc(2, 5, k)
            ti4 = cc(2, 3, k) - cc(2, 4, k)
            ti3 = cc(2, 3, k) + cc(2, 4, k)
            tr5 = cc(1, 2, k) - cc(1, 5, k)
            tr2 = cc(1, 2, k) + cc(1, 5, k)
            tr4 = cc(1, 3, k) - cc(1, 4, k)
            tr3 = cc(1, 3, k) + cc(1, 4, k)
            ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
            ch(2, k, 1) = cc(2, 1, k) + ti2 + ti3
            cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
            ci2 = cc(2, 1, k) + tr11*ti2 + tr12*ti3
            cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
            ci3 = cc(2, 1, k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            ch(1, k, 2) = cr2 - ci5
            ch(1, k, 5) = cr2 + ci5
            ch(2, k, 2) = ci2 + cr5
            ch(2, k, 3) = ci3 + cr4
            ch(1, k, 3) = cr3 - ci4
            ch(1, k, 4) = cr3 + ci4
            ch(2, k, 4) = ci3 - cr4
            ch(2, k, 5) = ci2 - cr5
         end do
      end if
   end subroutine passf5
end module fftpack_legacy_drivers_pass
