C
      subroutine orban(iorb,id,ar,br,n,l,occup,spin,eigv,ekin,epot)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
c      
c  orban is used to analyze and printout data
c  about the orbital.
c
C     .. Parameters ..
      double precision ai, zero
      parameter (ai=2*137.0360411D0,zero=0.D0)
C     ..
C     .. Scalar Arguments ..
      integer iorb, l, n
      double precision spin, occup, eigv, ekin, epot 
      character id*1
C     ..
C     .. Array Arguments ..
      double precision ar(*), br(*)
C     ..
C     .. Local Scalars ..
      double precision ar2, arp, arpm, br2, deni, sa2
      integer i, i90, i99, ka, ll, llp, lp, nextr, nzero
      integer kj, ist
C     ..
C     .. Local Arrays ..
      double precision aextr(10), bextr(10), rextr(10), rzero(10)
C     ..
      logical leqi
      external leqi
C     ..
c
      ka = l + 1
      lp = ka
      if (spin .lt. 0.1D0 .and. l .ne. 0) ka = -l
c
c      compute zeroes and extrema
c
      nzero = 0
      nextr = 0
      rzero(1) = zero
      arp = br(2)
      if (leqi(id,'r')) then
         if (spin .lt. 0.1D0) then
            arp = ka*ar(2)/r(2) + (eigv-viod(lp,2)/r(2)-vid(2)+
     &            ai*ai)*br(2)/ai
         else
            arp = ka*ar(2)/r(2) + (eigv-viou(lp,2)/r(2)-viu(2)+
     &            ai*ai)*br(2)/ai
         end if
      end if
      do 20 i = 3, nr
         if (nextr .ge. n-l) go to 30
c
         if (ar(i)*ar(i-1) .le. zero) then
c
c            zero
c
             nzero = nzero + 1
             rzero(nzero) = (ar(i)*r(i-1)-ar(i-1)*r(i))/
     &                                            (ar(i)-ar(i-1))
         endif
c
         arpm = arp
         arp = br(i)
c
         if (leqi(id,'r')) then
            if (spin .lt. 0.1D0) then
               arp = ka*ar(i)/r(i) + (eigv-viod(lp,i)/r(i)-vid(i)+
     &               ai*ai)*br(i)/ai
            else
               arp = ka*ar(i)/r(i) + (eigv-viou(lp,i)/r(i)-viu(i)+
     &               ai*ai)*br(i)/ai
            end if
         end if
c
         if (arp*arpm .le. zero) then
c
c          extremum
c
           nextr = nextr + 1
           rextr(nextr) = (arp*r(i-1)-arpm*r(i))/(arp-arpm)
           aextr(nextr) = (ar(i)+ar(i-1))/2 -
     &                    (arp**2+arpm**2)*(r(i)-r(i-1))/(4*(arp-arpm))
           bextr(nextr) = br(i)
c
         endif
c
   20 continue
c
c   Find orbital kinetic and potential energy
c   the potential part includes only the interaction with
c   the nuclear part
c
   30 continue
      ekin = br(1)*br(1)*rab(1)
      epot = zero
      sa2 = zero
      lp = l + 1
      llp = l*lp
      ll = 2
      if (2*(nr/2) .eq. nr) ll = 4
      i90 = nr
      i99 = nr
      do 40 i = nr, 2, -1
         ar2 = ar(i)*ar(i)
         br2 = br(i)*br(i)
         deni = ar2
         if (leqi(id,'r')) deni = deni + br2
         ekin = ekin + ll*(br2+ar2*llp/r(i)**2)*rab(i)
         if (spin .lt. 0.1D0) then
            epot = epot + ll*deni*viod(lp,i)*rab(i)/r(i)
         else
            epot = epot + ll*deni*viou(lp,i)*rab(i)/r(i)
         end if
         ll = 6 - ll
         if (sa2 .le. 0.1D0) then
           sa2 = sa2 + deni*rab(i)
           if (sa2 .le. 0.01D0) i99 = i
           i90 = i
         endif
   40 continue
c
      ekin = ekin/3
      epot = epot/3
      if (leqi(id,'r')) ekin = zero
c
c     Printout
c
      write(6,9000) n, l, spin
 9000 format(/' n =',i2,'  l =',i2,'  s =',f4.1)
c
      write(6,9010) 'a extr    ', (aextr(i),i=1,nextr)
      if (leqi(id,'r')) write(6,9010) 'b extr    ', (bextr(i),i=1,nextr)
      write(6,9010) 'r extr    ', (rextr(i),i=1,nextr)
      write(6,9010) 'r zero    ', (rzero(i),i=1,nzero)
      write(6,9010) 'r 90/99 % ', r(i90), r(i99)
c
      if (eigv .eq. zero) then
         if (occup .ne. zero) then
            write(6,9020) occup
         else
            write(6,9030)
         end if
      end if
c
 9010 format(8x,a10,2x,8f8.3)
 9020 format(8x,'WARNING: This orbital is not bound',' and contains ',
     &      f6.4,' electrons!!')
 9030 format(8x,'WARNING:  This orbital is not bound!')
c
c
      if (job .ne. 0 .and. job .ne. 4) return
c
c    Plot valence wavefunctions if AE or PT job
c
      if (iorb .gt. ncore) then
c
c     Plot and make the wavefunction 'upright'
c
         ist = nint(sign(1.d0,ar(nr-85)))
c
         if (job .eq. 4) then
            kj = -1
         else 
            kj = 1
         endif
         call potrw(ar,r,nr-85,l,kj,ist,0.d0)
      endif
c
      return
c
      end


