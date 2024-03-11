      subroutine cffti1(n,Wa,Ifac)
      use fftpack_kind
      implicit none
      real(rk) :: arg , argh , argld , fi , Wa
      integer :: i , i1 , ib , ido , idot , Ifac , ii , ip , ipm , j , k1, &
                 l1 , l2 , ld , n , nf , nl , nq , nr , ntry
      dimension Wa(*) , Ifac(*)  
      integer,dimension(4),parameter :: ntryh = [3 , 4 , 2 , 5]
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
      i = 2
      l1 = 1
      do k1 = 1 , nf
         ip = Ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido + ido + 2
         ipm = ip - 1
         do j = 1 , ipm
            i1 = i
            Wa(i-1) = 1.0_rk
            Wa(i) = 0.0_rk
            ld = ld + l1
            fi = 0.0_rk
            argld = real(ld, rk)*argh
            do ii = 4 , idot , 2
               i = i + 2
               fi = fi + 1.0_rk
               arg = fi*argld
               Wa(i-1) = cos(arg)
               Wa(i) = sin(arg)
            enddo
            if ( ip>5 ) then
               Wa(i1-1) = Wa(i-1)
               Wa(i1) = Wa(i)
            endif
         enddo
         l1 = l2
      enddo
      end subroutine cffti1