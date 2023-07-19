C
c $Id: vionic.f,v 1.4 2002/07/08 18:08:26 wdpgaara Exp $
c
c $Log: vionic.f,v $
c Revision 1.4  2002/07/08 18:08:26  wdpgaara
c Re-implemented the "valence charge modification" feature.
c New routine: change_valence.
c Superseded vionic kludge. Watch out for compiler complaints.
c
c Revision 1.3  2002/07/05 18:22:39  wdpgaara
c Fix format of grid parameters
c
c Revision 1.2  1997/05/22 17:32:36  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine Vionic
c
c  Vionic sets up the ionic potential.
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
c
c  njtj ***  major modifications  ***
c    If a potential does not exist, it is approximated
c    by an existing potential.
c    A nonspin or spin-polarized pseudo test, uses the
c    down(nonspin generation), weighted average(spin-
c    polarized), or averaged(relativistic) potentials.
c    A relativistic pseudo test, must use relativistic
c    generated potentials.  The Schroedinger equation is
c    used to integrate a relativistic pseudo test,
c    not the Dirac equation.
c  njtj  ***  major modifications  ***
c
C     .. Parameters ..
      double precision zero
      parameter (zero=0.D0)
C     ..
C     .. Local Scalars ..
      double precision vdiff, vsum, zion
      integer i, j, loi, npotd, npotu, nrm
      character icorrt*2, namet*2
C     ..
C     .. Local Arrays ..
      integer npd(5), npu(5)
      character ray(6)*10, title(7)*10
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
c
c  2*znuc part
c
      ifcore = 0
      if (job .lt. 4) then
         do 20 i = 1, lmax
            do 10 j = 1, nr
               viod(i,j) = -2*znuc
               viou(i,j) = -2*znuc
   10       continue
   20    continue
      else
c
c  read pseudopotentials from tape1
c
         open(unit=1,file='VPSIN',status='old',form='unformatted')
         rewind 1
         read(1) namet, icorrt, irel, nicore, (ray(i),i=1,6),
     &     (title(i),i=1,7), npotd, npotu, nrm, a, b, zion
c
         if (nicore .eq. 'fcec' .or. nicore .eq. 'pcec') ifcore = 1
         if (nicore .eq. 'fche' .or. nicore .eq. 'pche') ifcore = 2
         nr = nrm + 1
         read(1) (r(i),i=2,nr)
         r(1) = zero
c
c   down potentials (or average relativistic potentials)
c
c njtj  ***  major start  ***
c   if a potential does not exist, it is replaced by the
c   next existing lower angular momentum potential or
c   the next existing higher if no lower exists.
c
         do 30 i = 1, lmax
            npd(i) = 0
   30    continue
         do 40 i = 1, npotd
            read(1) loi, (viod(loi+1,j),j=2,nr)
            viod(loi+1,1) = zero
            npd(loi+1) = 1
   40    continue
         if (npd(1) .eq. 0) then
            do 60 i = 2, lmax
               if (npd(i) .gt. 0) then
                  do 50 j = 1, nr
                     viod(1,j) = viod(i,j)
   50             continue
c
                  go to 70
c
               end if
   60       continue
         end if
   70    continue
         do 90 i = 2, lmax
            if (npd(i) .eq. 0) then
               do 80 j = 1, nr
                  viod(i,j) = viod(i-1,j)
   80          continue
            end if
   90    continue
c
c   up potentials (or spin orbit potentials)
c
         if (npotu .le. 0) go to 170
         do 100 i = 1, lmax
            npu(i) = 0
  100    continue
         do 110 i = 1, npotu
            read(1) loi, (viou(loi+1,j),j=2,nr)
            viou(loi+1,1) = zero
            npu(loi+1) = 1
  110    continue
         if (npu(1) .eq. 0) then
            do 130 i = 2, lmax
               if (npu(i) .gt. 0) then
                  do 120 j = 1, nr
                     viou(1,j) = viou(i,j)
  120             continue
c
                  go to 140
c
               end if
  130       continue
         end if
  140    continue
         do 160 i = 2, lmax
            if (npu(i) .eq. 0) then
               do 150 j = 1, nr
                  viou(i,j) = viou(i-1,j)
  150          continue
            end if
  160    continue
c
c  njtj  ***  major end  ***
c
c
c  core and valence charges
c
  170    continue
         read(1) (cdc(i),i=2,nr)
         cdc(1) = zero
c
c  replace valence charge on tape(valence charge modify)
c  This is a horrible kludge, and it probably has never
c  been used after the Froyen days...
c  It might also violate the Fortran Standard.
c
c  Note the instantaneous job number change needed
c  Leave it as is, but use routine change_valence
c  for a cleaner operation.
c
         if (job .eq. 6) then
            write(1) (cdd(i)+cdu(i),i=2,nr)
            close(1)
c
            return
c
         end if
         read(1) (cdd(i),i=2,nr)

         close(1)

         cdd(1) = zero
c
c  njtj  ***   major start  ***
c   distribute charge as up and down charge
c   generate radial intergration grid
c   set up potentials equal to down potentials for
c   spin-polarized pseudo test of nonspin and relativistic
c   generated potentails.  Construct spin-orbit potentials
c   from relativistic sum and difference potentials and
c   change ispp='r' to ispp=' '.
c
         do 180 i = 1, nr
            rab(i) = (r(i)+a)*b
            cdd(i) = cdd(i)/2
            cdu(i) = cdd(i)
  180    continue
         if (ispp .eq. 's' .and. irel .ne. 'isp') then
            do 200 i = 1, lmax
               do 190 j = 1, nr
                  viou(i,j) = viod(i,j)
  190          continue
  200       continue
         end if
         if (ispp .eq. 'r') then
            ispp = ' '
            if (irel .ne. 'rel') then
               write(6,9000)
 9000          format(//'Pseudopotential is not relativistic!!!!',
     &               /' setting up potentials equal to down!!!',//)
               do 220 i = 1, lmax
                  do 210 j = 1, nr
                     viou(i,j) = viod(i,j)
  210             continue
  220          continue
            else
               do 230 j = 1, nr
                  viou(1,j) = viod(1,j)
  230          continue
               do 250 i = 2, lmax
                  do 240 j = 1, nr
                     vsum = viod(i,j)
                     vdiff = viou(i,j)
                     viod(i,j) = vsum - i*vdiff/2
                     viou(i,j) = vsum + (i-1)*vdiff/2
  240             continue
  250          continue
            end if
         end if
c
c   njtj  ***  major end   ***
c
c
c   printout
c
         write(6,9010) namet, icorrt, irel, nicore,
     &     (ray(i),i=1,6), (title(i),i=1,7)
 9010    format(//1x,a2,2x,a2,2x,a3,2x,a4,
     &         '  pseudopotential read from tape',/1x,2a10,5x,4a10,/1x,
     &         7a10,//)
         if (nameat .ne. namet) write(6,9020) nameat, namet
 9020    format(' input element ',a2,' not equal to element on tape ',
     &         a2,//)
         if (icorr .ne. icorrt) write(6,9030) icorr, icorrt
 9030    format(' input correlation ',a2,
     &         ' not equal to correlation from tape ',a2,//)
         write(6,9040) r(2), nr, r(nr)
 9040    format(' radial grid parameters',//' r(1) = .0 , r(2) =',d8.2,
     &         ' , ... , r(',i4,') =',f6.2,//)
      end if
c
c   add potential from shell charge
c
      if (abs(zsh) .gt. 0.D-5) then
         do 270 i = 1, lmax
            do 260 j = 1, nr
               if (r(j) .ge. rsh) then
                  viod(i,j) = viod(i,j) - 2*zsh
                  viou(i,j) = viou(i,j) - 2*zsh
               else
                  viod(i,j) = viod(i,j) - 2*zsh*r(j)/rsh
                  viou(i,j) = viou(i,j) - 2*zsh*r(j)/rsh
               end if
  260       continue
  270    continue
      end if
c
      return
c
      end
