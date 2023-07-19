C
      subroutine change_valence
c
c     Generates a new pseudopotential file
c     with a modified valence charge in it.
c     This might be useful ...
c
      include 'radial.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
c
C     .. Local Scalars ..
      double precision zion
      integer i, j, loi, npotd, npotu, nrm
      character icorrt*2, namet*2
C     ..
C     .. Local Arrays ..
      character ray(6)*10, title(7)*10
C     ..


      open(unit=1,file='VPSIN',status='old',form='unformatted')
      open(unit=2,file='VPS_NEW_VALENCE',
     $        status='unknown',form='unformatted')
      rewind 1
      rewind 2
      read(1) namet, icorrt, irel, nicore, (ray(i),i=1,6),
     &     (title(i),i=1,7), npotd, npotu, nrm, a, b, zion
c
      read(1) (r(i),i=2,nr)
      ray(6) = 'NEWVALENCE'
      write(2) namet, icorrt, irel, nicore, (ray(i),i=1,6),
     &     (title(i),i=1,7), npotd, npotu, nrm, a, b, zion
c
      write(2) (r(i),i=2,nr)
c
c   down potentials (or average relativistic potentials)
c
      do 40 i = 1, npotd
         read(1) loi, (viod(loi+1,j),j=2,nr)
         write(2) loi, (viod(loi+1,j),j=2,nr)
 40   continue
c
c   up potentials (or spin orbit potentials)
c
      if (npotu .gt. 0) then
         do 110 i = 1, npotu
            read(1) loi, (viou(loi+1,j),j=2,nr)
            write(2) loi, (viou(loi+1,j),j=2,nr)
  110    continue
      endif
c     
c  core and valence charges
c
      read(1) (cdc(i),i=2,nr)
      write(2) (cdc(i),i=2,nr)

c     Write the actual valence charge we
c     have in the program's structures
c     Without the infamous "zratio" ...
c
      write(2) (cdd(i)+cdu(i),i=2,nr)

      close(1)
      close(2)

      end






