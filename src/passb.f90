!*==PASSB.spg  processed by SPAG 6.72Dc at 19:17 on 14 Sep 2021
      subroutine passb(Nac,Ido,Ip,l1,Idl1,Cc,c1,c2,Ch,Ch2,Wa)
      use fftpack_kind
      implicit none
!*--PASSB722
!*** Start of declarations inserted by SPAG
      real c1 , c2 , Cc , Ch , Ch2 , fftpack_kind , rk , Wa , wai , war
      integer i , idij , idj , idl , Idl1 , idlj , Ido , idot , idp ,   &
            & ik , inc , Ip , ipp2 , ipph , j , jc , k , l , l1 , lc
      integer Nac , nt
!*** End of declarations inserted by SPAG
      dimension Ch(Ido,l1,Ip) , Cc(Ido,Ip,l1) , c1(Ido,l1,Ip) , Wa(1) , &
              & c2(Idl1,Ip) , Ch2(Idl1,Ip)
      idot = Ido/2
      nt = Ip*Idl1
      ipp2 = Ip + 2
      ipph = (Ip+1)/2
      idp = Ip*Ido
!
      if ( Ido<l1 ) then
         do j = 2 , ipph
            jc = ipp2 - j
            do i = 1 , Ido
               do k = 1 , l1
                  Ch(i,k,j) = Cc(i,j,k) + Cc(i,jc,k)
                  Ch(i,k,jc) = Cc(i,j,k) - Cc(i,jc,k)
               enddo
            enddo
         enddo
         do i = 1 , Ido
            do k = 1 , l1
               Ch(i,k,1) = Cc(i,1,k)
            enddo
         enddo
      else
         do j = 2 , ipph
            jc = ipp2 - j
            do k = 1 , l1
               do i = 1 , Ido
                  Ch(i,k,j) = Cc(i,j,k) + Cc(i,jc,k)
                  Ch(i,k,jc) = Cc(i,j,k) - Cc(i,jc,k)
               enddo
            enddo
         enddo
         do k = 1 , l1
            do i = 1 , Ido
               Ch(i,k,1) = Cc(i,1,k)
            enddo
         enddo
      endif
      idl = 2 - Ido
      inc = 0
      do l = 2 , ipph
         lc = ipp2 - l
         idl = idl + Ido
         do ik = 1 , Idl1
            c2(ik,l) = Ch2(ik,1) + Wa(idl-1)*Ch2(ik,2)
            c2(ik,lc) = Wa(idl)*Ch2(ik,Ip)
         enddo
         idlj = idl
         inc = inc + Ido
         do j = 3 , ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if ( idlj>idp ) idlj = idlj - idp
            war = Wa(idlj-1)
            wai = Wa(idlj)
            do ik = 1 , Idl1
               c2(ik,l) = c2(ik,l) + war*Ch2(ik,j)
               c2(ik,lc) = c2(ik,lc) + wai*Ch2(ik,jc)
            enddo
         enddo
      enddo
      do j = 2 , ipph
         do ik = 1 , Idl1
            Ch2(ik,1) = Ch2(ik,1) + Ch2(ik,j)
         enddo
      enddo
      do j = 2 , ipph
         jc = ipp2 - j
         do ik = 2 , Idl1 , 2
            Ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
            Ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
            Ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
            Ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
         enddo
      enddo
      Nac = 1
      if ( Ido==2 ) return
      Nac = 0
      do ik = 1 , Idl1
         c2(ik,1) = Ch2(ik,1)
      enddo
      do j = 2 , Ip
         do k = 1 , l1
            c1(1,k,j) = Ch(1,k,j)
            c1(2,k,j) = Ch(2,k,j)
         enddo
      enddo
      if ( idot>l1 ) then
         idj = 2 - Ido
         do j = 2 , Ip
            idj = idj + Ido
            do k = 1 , l1
               idij = idj
               do i = 4 , Ido , 2
                  idij = idij + 2
                  c1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) - Wa(idij)       &
                              & *Ch(i,k,j)
                  c1(i,k,j) = Wa(idij-1)*Ch(i,k,j) + Wa(idij)           &
                            & *Ch(i-1,k,j)
               enddo
            enddo
         enddo
         goto 99999
      endif
      idij = 0
      do j = 2 , Ip
         idij = idij + 2
         do i = 4 , Ido , 2
            idij = idij + 2
            do k = 1 , l1
               c1(i-1,k,j) = Wa(idij-1)*Ch(i-1,k,j) - Wa(idij)*Ch(i,k,j)
               c1(i,k,j) = Wa(idij-1)*Ch(i,k,j) + Wa(idij)*Ch(i-1,k,j)
            enddo
         enddo
      enddo
      return
99999 end subroutine passb