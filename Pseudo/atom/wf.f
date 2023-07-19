c
      subroutine wf(i,ar,br)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
c     Solves the wave equation for orbital i and
c     finds the extrema and nodes.
c
      double precision zero, ai, pnine
      parameter (zero=0.d0,ai=2*137.0360411D0,pnine=0.9d0)
c
      integer i, nextr
      double precision rextr, rzero
c
      double precision ar(*), br(*)
c
      integer iflag, ist, j, lp, llp, ka
      double precision arp, arpm
      double precision v(nrmax)
c
      lp = lo(i) + 1
      llp = lo(i)*lp
c
      do 10 j = 1, nr
         ar(j) = 0.d0
   10 continue
      if (down(i)) then
         do 20 j = 2, nr
            v(j) = viod(lp,j)/r(j) + vid(j)
   20    continue
      else
         do 30 j = 2, nr
            v(j) = viou(lp,j)/r(j) + viu(j)
   30    continue
      end if
c
      if ( .not. relativistic) then
c
c           Add 'centrifugal term'
c
         do 40 j = 2, nr
            v(j) = v(j) + llp/r(j)**2
   40    continue
      end if
c
      if ( .not. relativistic) then
         call difnrl(0,i,v,ar,br,no(i),lo(i),so(i),ev(i),iflag)
      else
         call difrel(0,i,v,ar,br,no(i),lo(i),so(i),ev(i))
      end if
c
c     Plot and make the wavefunction 'upright'
c
      ist = nint(sign(1.d0,ar(nr-85)))
      call potrw(ar,r,nr-85,lo(i),1,ist,rc(lp))
c
      do 50 j = 1, nr
         ar(j) = ar(j)*ist
         br(j) = br(j)*ist
   50 continue
c
c  Find the last zero and extremum.
c
      ka = lo(i) + 1
      if (down(i) .and. lo(i) .ne. 0) ka = -lo(i)
      nextr = no(i) - lo(i)
      rzero = zero
      arp = br(2)
c
      if (relativistic) then
         if (down(i)) then
            arp = ka*ar(2)/r(2) + (ev(i)-viod(lp,2)/r(2)-vid(2)+ai*ai)*
     &            br(2)/ai
         else
            arp = ka*ar(2)/r(2) + (ev(i)-viou(lp,2)/r(2)-viu(2)+ai*ai)*
     &            br(2)/ai
         end if
      end if
c
      do 60 j = 3, nr
c
         if (nextr .eq. 0) go to 70
c
         if (ar(j-1)*ar(j) .le. zero) rzero = (ar(j)*r(j-1)-
     &       ar(j-1)*r(j))/(ar(j)-ar(j-1))
         arpm = arp
         arp = br(j)
c
         if (relativistic) then
            if (down(i)) then
               arp = ka*ar(j)/r(j) + (ev(i)-viod(lp,j)/r(j)-vid(j)+
     &               ai*ai)*br(j)/ai
            else
               arp = ka*ar(j)/r(j) + (ev(i)-viou(lp,j)/r(j)-viu(j)+
     &               ai*ai)*br(j)/ai
            end if
         end if
c
         if (arp*arpm .le. zero) then
            rextr = (arp*r(j-1)-arpm*r(j))/(arp-arpm)
            nextr = nextr - 1
         end if
c
   60 continue
   70 continue
c
c  Check rc, if outside bounds reset.
c
      if (rzero .lt. r(2)) rzero = r(2)
c
c  Check rc if inside rzero,
c  reset to .9 between rmax and rzero if inside
c  if rc(lp) is negative, rc(lp) is percent of way
c  betweeen rzero and rmax.
c
      if (rc_input(lp) .gt. rzero) then
c
c           do nothing
c
      else if (rc_input(lp) .ge. zero) then
c
c        rc is inside the node...
c
         write(6,'(a,3f8.3)')
     $     ' Requested rc inside node ! ** rzero, rextr, rc_input:',
     $                          rzero, rextr, rc_input(lp)
         rc_input(lp) = rzero + pnine*(rextr-rzero)
         write(6,'(a,f6.3)') ' rc changed to ', rc_input(lp)
c
      else
c
c        compute rc as a fraction of rextr-rzero
c
         rc_input(lp) = rzero - rc_input(lp)*(rextr-rzero)
         write(6,'(a,f5.2)') 'rc set to ', rc_input(lp)
c
      end if
c
c        Warn the user if rc > rextr
c
      if (rc_input(lp) .ge. rextr) write(6,9000) rc_input(lp), rextr
 9000 format(' Core radius (',f5.2,
     &      ') outside wfn extremum (',f5.2,')')
c
      return
c
      end


