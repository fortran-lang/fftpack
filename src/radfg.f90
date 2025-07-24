      subroutine radfg(Ido,Ip,l1,Idl1,Cc,c1,c2,Ch,Ch2,Wa)
      use fftpack_kind
      implicit none
      real(rk) :: ai1 , ai2 , ar1 , ar1h , ar2 , ar2h , arg , c1 , &
                  c2 , Cc , Ch , Ch2 , dc2 , dcp , ds2 , dsp , &
                  Wa
      integer :: i , ic , idij , Idl1 , Ido , idp2 , ik , Ip , ipp2 , &
                 ipph , is , j , j2 , jc , k , l , l1 , lc , nbd
      dimension Ch(Ido,l1,Ip) , Cc(Ido,Ip,l1) , c1(Ido,l1,Ip) , &
                c2(Idl1,Ip) , Ch2(Idl1,Ip) , Wa(*)
      real(rk),parameter :: tpi = 2.0_rk * acos(-1.0_rk) ! 2 * pi
      arg = tpi/real(Ip, rk)
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (Ip+1)/2
      ipp2 = Ip + 2
      idp2 = Ido + 2
      nbd = (Ido-1)/2
      if ( Ido==1 ) then
         do ik = 1 , Idl1
            c2(ik,1) = Ch2(ik,1)
         enddo
      else
         do ik = 1 , Idl1
            Ch2(ik,1) = c2(ik,1)
         enddo
         do j = 2 , Ip
            do k = 1 , l1
               Ch(1,k,j) = c1(1,k,j)
            enddo
         enddo
         if ( nbd>l1 ) then
            is = -Ido
            do j = 2 , Ip
               is = is + Ido
               do k = 1 , l1
                  idij = is
                  do i = 3 , Ido , 2
                     idij = idij + 2
                     Ch(i-1,k,j) = Wa(idij-1)*c1(i-1,k,j) + Wa(idij) &
                                   *c1(i,k,j)
                     Ch(i,k,j) = Wa(idij-1)*c1(i,k,j) - Wa(idij) &
                                 *c1(i-1,k,j)
                  enddo
               enddo
            enddo
         else
            is = -Ido
            do j = 2 , Ip
               is = is + Ido
               idij = is
               do i = 3 , Ido , 2
                  idij = idij + 2
                  do k = 1 , l1
                     Ch(i-1,k,j) = Wa(idij-1)*c1(i-1,k,j) + Wa(idij) &
                                   *c1(i,k,j)
                     Ch(i,k,j) = Wa(idij-1)*c1(i,k,j) - Wa(idij) &
                                 *c1(i-1,k,j)
                  enddo
               enddo
            enddo
         endif
         if ( nbd<l1 ) then
            do j = 2 , ipph
               jc = ipp2 - j
               do i = 3 , Ido , 2
                  do k = 1 , l1
                     c1(i-1,k,j) = Ch(i-1,k,j) + Ch(i-1,k,jc)
                     c1(i-1,k,jc) = Ch(i,k,j) - Ch(i,k,jc)
                     c1(i,k,j) = Ch(i,k,j) + Ch(i,k,jc)
                     c1(i,k,jc) = Ch(i-1,k,jc) - Ch(i-1,k,j)
                  enddo
               enddo
            enddo
         else
            do j = 2 , ipph
               jc = ipp2 - j
               do k = 1 , l1
                  do i = 3 , Ido , 2
                     c1(i-1,k,j) = Ch(i-1,k,j) + Ch(i-1,k,jc)
                     c1(i-1,k,jc) = Ch(i,k,j) - Ch(i,k,jc)
                     c1(i,k,j) = Ch(i,k,j) + Ch(i,k,jc)
                     c1(i,k,jc) = Ch(i-1,k,jc) - Ch(i-1,k,j)
                  enddo
               enddo
            enddo
         endif
      endif
      do j = 2 , ipph
         jc = ipp2 - j
         do k = 1 , l1
            c1(1,k,j) = Ch(1,k,j) + Ch(1,k,jc)
            c1(1,k,jc) = Ch(1,k,jc) - Ch(1,k,j)
         enddo
      enddo
!
      ar1 = 1.0_rk
      ai1 = 0.0_rk
      do l = 2 , ipph
         lc = ipp2 - l
         ar1h = dcp*ar1 - dsp*ai1
         ai1 = dcp*ai1 + dsp*ar1
         ar1 = ar1h
         do ik = 1 , Idl1
            Ch2(ik,l) = c2(ik,1) + ar1*c2(ik,2)
            Ch2(ik,lc) = ai1*c2(ik,Ip)
         enddo
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do j = 3 , ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h
            do ik = 1 , Idl1
               Ch2(ik,l) = Ch2(ik,l) + ar2*c2(ik,j)
               Ch2(ik,lc) = Ch2(ik,lc) + ai2*c2(ik,jc)
            enddo
         enddo
      enddo
      do j = 2 , ipph
         do ik = 1 , Idl1
            Ch2(ik,1) = Ch2(ik,1) + c2(ik,j)
         enddo
      enddo
!
      if ( Ido<l1 ) then
         do i = 1 , Ido
            do k = 1 , l1
               Cc(i,1,k) = Ch(i,k,1)
            enddo
         enddo
      else
         do k = 1 , l1
            do i = 1 , Ido
               Cc(i,1,k) = Ch(i,k,1)
            enddo
         enddo
      endif
      do j = 2 , ipph
         jc = ipp2 - j
         j2 = j + j
         do k = 1 , l1
            Cc(Ido,j2-2,k) = Ch(1,k,j)
            Cc(1,j2-1,k) = Ch(1,k,jc)
         enddo
      enddo
      if ( Ido==1 ) return
      if ( nbd<l1 ) then
         do j = 2 , ipph
            jc = ipp2 - j
            j2 = j + j
            do i = 3 , Ido , 2
               ic = idp2 - i
               do k = 1 , l1
                  Cc(i-1,j2-1,k) = Ch(i-1,k,j) + Ch(i-1,k,jc)
                  Cc(ic-1,j2-2,k) = Ch(i-1,k,j) - Ch(i-1,k,jc)
                  Cc(i,j2-1,k) = Ch(i,k,j) + Ch(i,k,jc)
                  Cc(ic,j2-2,k) = Ch(i,k,jc) - Ch(i,k,j)
               enddo
            enddo
         enddo
      else
         do j = 2 , ipph
            jc = ipp2 - j
            j2 = j + j
            do k = 1 , l1
               do i = 3 , Ido , 2
                  ic = idp2 - i
                  Cc(i-1,j2-1,k) = Ch(i-1,k,j) + Ch(i-1,k,jc)
                  Cc(ic-1,j2-2,k) = Ch(i-1,k,j) - Ch(i-1,k,jc)
                  Cc(i,j2-1,k) = Ch(i,k,j) + Ch(i,k,jc)
                  Cc(ic,j2-2,k) = Ch(i,k,jc) - Ch(i,k,j)
               enddo
            enddo
         enddo
      end if
      end subroutine radfg