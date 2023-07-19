c
      subroutine dsolv2(iter,iconv,id,nfirst,nlast,n_of_core_orbs,nn)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
c
c  dsolv2 finds the (non) relativistic wave function using
c  difnrl to integrate the Schroedinger equation or
c  difrel to integrate the Dirac equation.
c  The energy level from the previous iteration are used
c  as initial guesses, and it must therefore be reasonably
c  accurate.
c
C     .. Parameters ..
      double precision zero, smev
      parameter (zero=0.D0,smev=1.D-4)
C     ..
C     .. Scalar Arguments ..
      integer iconv, iter, nfirst, nlast, n_of_core_orbs
      character id*1
C     ..
c
c     Watch out for nn: it is not always equal to 'no'. 
c     In particular, some nn(i) could be zero.
c     The same is true for id...
c
      integer nn(*)
C     ..
C     .. Local Scalars ..
      integer i, iflag, j, llp, lp
C     ..
C     .. External Subroutines ..
      external difnrl, difrel, orban
      logical leqi
      external leqi
C     ..
      double precision ar(nrmax), br(nrmax), v(nrmax),
     &                 orb_charge(nrmax)
C     ..
c
c  Initialize arrays for charge density.
c
      do 10 i = 1, nr
         cdd(i) = zero
         cdu(i) = zero
   10 continue
      if (ifcore .ne. 1) then
         do 30 i = 1, nr
            cdc(i) = zero
   30    continue
      end if
c
c  Start the loop over orbitals.
c  Note that spin zero is treated as down.
c
      do 130 i = nfirst, nlast
         if (nn(i) .le. 0) go to 130
         if (zo(i) .eq. 0.0D0 .and. iconv .eq. 0) go to 130
         if (ev(i) .ge. 0.0D0) ev(i) = -smev
c
c  Set up the potential, set the wave function array to zero-ar.
c
         lp = lo(i) + 1
         llp = lo(i)*lp
         do 40 j = 1, nr
            ar(j) = zero
   40    continue
         if (down(i)) then
            do 50 j = 2, nr
               v(j) = viod(lp,j)/r(j) + vid(j)
   50       continue
         else
            do 60 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j)
   60       continue
         end if
         if (.not. leqi(id,'r')) then
            do 70 j = 2, nr
               v(j) = v(j) + llp/r(j)**2
   70       continue
         end if
c
c  Call the integration routine.
c
         if (.not. leqi(id,'r')) then
            call difnrl(iter,i,v,ar,br,nn(i),lo(i),so(i),ev(i),iflag)
         else
            call difrel(iter,i,v,ar,br,nn(i),lo(i),so(i),ev(i))
         end if
c
c  Add to the charge density.
c
         if (leqi(id,'r')) then
            do 300 j = 1, nr
               orb_charge(j) = zo(i)*(br(j)*br(j)+ar(j)*ar(j))
 300        continue
         else
            do 320 j = 1, nr
               orb_charge(j) = zo(i) * ar(j)*ar(j)
 320        continue
         endif
c
            if (down(i)) then
               do 80 j = 1, nr
                  cdd(j) = cdd(j) + orb_charge(j)
   80          continue
            else
               do 90 j = 1, nr
                  cdu(j) = cdu(j) + orb_charge(j)
   90          continue
            end if

cag            if (ifcore .ne. 1 .and. i .le. n_of_core_orbs) then
            if ( i .le. n_of_core_orbs) then
              do 95 j = 1, nr
                 cdc(j) = cdc(j) + orb_charge(j)
   95         continue
            end if

c
c  Compute various quantitities if last iteration.
c
         if (iconv .eq. 1) call orban(i,id,ar,br,nn(i),lo(i),zo(i),
     &                                so(i),ev(i),ek(i),ep(i))
  130 continue
c
c  End loop over orbitals.
c
      return
c
      end
