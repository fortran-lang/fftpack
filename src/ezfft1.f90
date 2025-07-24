      subroutine ezfft1(n,Wa,Ifac)
      use fftpack_kind
      implicit none
      real(rk) :: arg1 , argh , ch1 , ch1h , dch1 , dsh1 , sh1 , &
                  Wa
      integer :: i , ib , ido , Ifac , ii , ip , ipm , is , j , k1 , l1 , &
                 l2 , n , nf , nfm1 , nl , nq , nr , ntry
      dimension Wa(*) , Ifac(*)  
      integer,dimension(4),parameter :: ntryh = [4 , 2 , 3 , 5]
      real(rk),parameter :: tpi = 2.0_rk * acos(-1.0_rk) ! 2 * pi
      nl = n
      nf = 0
      j = 0
 100  j = j + 1
      if ( j<=4 ) then
         ntry = ntryh(j)
      else
         ntry = ntry + 2
      endif
 200  nq = nl/ntry
      nr = nl - ntry*nq
      if ( nr/=0 ) goto 100
      nf = nf + 1
      Ifac(nf+2) = ntry
      nl = nq
      if ( ntry==2 ) then
         if ( nf/=1 ) then
            do i = 2 , nf
               ib = nf - i + 2
               Ifac(ib+2) = Ifac(ib+1)
            enddo
            Ifac(3) = 2
         endif
      endif
      if ( nl/=1 ) goto 200
      Ifac(1) = n
      Ifac(2) = nf
      argh = tpi/real(n, rk)
      is = 0
      nfm1 = nf - 1
      l1 = 1
      if ( nfm1==0 ) return
      do k1 = 1 , nfm1
         ip = Ifac(k1+2)
         l2 = l1*ip
         ido = n/l2
         ipm = ip - 1
         arg1 = real(l1, rk)*argh
         ch1 = 1.0_rk
         sh1 = 0.0_rk
         dch1 = cos(arg1)
         dsh1 = sin(arg1)
         do j = 1 , ipm
            ch1h = dch1*ch1 - dsh1*sh1
            sh1 = dch1*sh1 + dsh1*ch1
            ch1 = ch1h
            i = is + 2
            Wa(i-1) = ch1
            Wa(i) = sh1
            if ( ido>=5 ) then
               do ii = 5 , ido , 2
                  i = i + 2
                  Wa(i-1) = ch1*Wa(i-3) - sh1*Wa(i-2)
                  Wa(i) = ch1*Wa(i-2) + sh1*Wa(i-3)
               enddo
            endif
            is = is + ido
         enddo
         l1 = l2
      enddo
      end subroutine ezfft1